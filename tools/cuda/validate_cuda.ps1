<#
.SYNOPSIS
    End-to-end CUDA validation: fresh CPU-only and CUDA Release builds, the
    full CTest suite, the CUDA-labelled tests, direct runs of the CUDA test
    binaries, Compute Sanitizer and a repeated-run stability check.

.DESCRIPTION
    Each stage is reported separately in <OutDir>/validation-summary.md and
    validation-summary.json, with full logs next to them.

    Test outcomes are classified, not just counted:
      Passed / Skipped      as CTest reports them
      Failed                the test ran and failed
      KnownFailure          failed, and listed in -KnownFailures (pre-existing,
                            unrelated issues; still reported)
      LaunchBlocked         CTest could not start the process (BAD_COMMAND /
                            Not Run). Every such test is retried by running its
                            command directly, and the direct result is
                            recorded: an operating-system launch block (for
                            example Windows Application Control) is reported
                            as such, never as a pass and never as a solver
                            failure.

    PYTHONUTF8=1 is set for test runs so the Python-driven tests decode child
    output as UTF-8 instead of the Windows ANSI code page.

    Exit code: 0 when every stage passed (known failures and OS launch blocks
    are reported but do not fail the run), 1 otherwise.

.EXAMPLE
    # Everything, on a machine with an sm_120 GPU
    tools/cuda/validate_cuda.ps1 -CudaArchitectures 120

.EXAMPLE
    # Reuse existing build trees, skip sanitizers
    tools/cuda/validate_cuda.ps1 -SkipBuild -SkipSanitizers
#>
param(
    [string] $CpuBuildDir = 'build-validate-cpu',
    [string] $CudaBuildDir = 'build-validate-cuda',
    [string] $CudaArchitectures = 'native',
    [string] $OutDir = 'build-validate-logs',
    [int] $DeterminismRuns = 10,
    [string[]] $KnownFailures = @('nlp_elastic_kkt_tests'),
    [switch] $SkipBuild,
    [switch] $SkipCpuBuild,
    [switch] $SkipCpuTests,
    [switch] $SkipSanitizers,
    [switch] $SkipDeterminism,
    [switch] $Racecheck
)

$ErrorActionPreference = 'Stop'
Import-Module (Join-Path $PSScriptRoot 'CudaValidation.psm1') -Force
$root = Get-RepoRoot
Push-Location $root
try {

$OutDir = New-OutputDirectory $OutDir
$testEnvironment = @{ PYTHONUTF8 = '1' }
$stages = [ordered]@{}
$overallOk = $true

function Add-Stage([string] $Name, [string] $Result, [string] $Detail, [bool] $Ok) {
    $stages[$Name] = [pscustomobject]@{ Result = $Result; Detail = $Detail; Ok = $Ok }
    if (-not $Ok) { $script:overallOk = $false }
    Write-Host ("[{0}] {1}: {2}" -f $(if ($Ok) { ' ok ' } else { 'FAIL' }), $Name, $Result)
    if ($Detail) { Write-Host "       $($Detail -replace "`n", "`n       ")" }
}

# --- CTest with classification and direct retries ---------------------------
function Invoke-CTestStage {
    param([string] $Name, [string] $BuildDir, [string[]] $ExtraArgs = @())
    $ctest = (Get-Command ctest -ErrorAction Stop).Source
    $logDir = New-OutputDirectory (Join-Path $OutDir $Name)
    $run = Invoke-Native -FilePath $ctest -Environment $testEnvironment -LogPath (Join-Path $logDir 'ctest.log') `
        -Arguments (@('--test-dir', $BuildDir, '-C', 'Release', '--output-on-failure', '--timeout', '1800') + $ExtraArgs)
    if (-not $run.Launched) {
        Add-Stage $Name 'ctest could not start' $run.LaunchError $false
        return
    }
    $tests = @()
    foreach ($line in ($run.Output -split "`r?`n")) {
        if ($line -match 'Test\s+#\d+:\s+(\S+)\s+\.*\s*(\*{3})?(.*?)\s+[0-9.]+\s+sec') {
            # Copy the captures first: switch -Regex below overwrites $Matches.
            $testName = $Matches[1]
            $status = $Matches[3].Trim()
            $category = switch -Regex ($status) {
                '^Passed' { 'Passed'; break }
                'Skipped' { 'Skipped'; break }
                'Not Run|BAD_COMMAND' { 'LaunchBlocked'; break }
                'Timeout' { 'Timeout'; break }
                default { 'Failed' }
            }
            if ($category -eq 'Failed' -and $KnownFailures -contains $testName) { $category = 'KnownFailure' }
            $tests += [pscustomobject]@{ Name = $testName; Status = $status; Category = $category; Direct = $null }
        }
    }

    # Retry every test CTest could not launch by running its command directly.
    $blocked = @($tests | Where-Object { $_.Category -eq 'LaunchBlocked' })
    if ($blocked.Count -gt 0) {
        $show = Invoke-Native -FilePath $ctest -Arguments @('--test-dir', $BuildDir, '-C', 'Release', '--show-only=json-v1')
        $definitions = ($show.Output | ConvertFrom-Json).tests
        foreach ($test in $blocked) {
            $definition = $definitions | Where-Object { $_.name -eq $test.Name } | Select-Object -First 1
            if (-not $definition -or -not $definition.command) { $test.Direct = 'no command recorded'; continue }
            $workDir = $BuildDir
            $env = @{} + $testEnvironment
            foreach ($property in @($definition.properties)) {
                if ($property.name -eq 'WORKING_DIRECTORY') { $workDir = $property.value }
                if ($property.name -eq 'ENVIRONMENT') {
                    foreach ($pair in @($property.value)) { $kv = $pair -split '=', 2; $env[$kv[0]] = $kv[1] }
                }
            }
            $command = @($definition.command)
            $direct = Invoke-Native -FilePath $command[0] -Arguments @($command | Select-Object -Skip 1) `
                -WorkingDirectory $workDir -Environment $env -LogPath (Join-Path $logDir "$($test.Name).direct.log")
            $verdict = Get-RunVerdict $direct @(77)
            $test.Direct = if ($verdict -eq 'LaunchBlocked') { "LaunchBlocked: $($direct.LaunchError)" } else { "$verdict (exit $($direct.ExitCode))" }
            if ($verdict -eq 'Passed') { $test.Category = 'PassedDirect' }
            elseif ($verdict -eq 'Skipped') { $test.Category = 'Skipped' }
            elseif ($verdict -eq 'Failed') {
                $test.Category = if ($KnownFailures -contains $test.Name) { 'KnownFailure' } else { 'Failed' }
            }
        }
    }

    $tests | ConvertTo-Json -Depth 3 | Set-Content -Path (Join-Path $logDir 'tests.json') -Encoding UTF8
    $count = { param($c) @($tests | Where-Object { $_.Category -eq $c }).Count }
    $summary = "{0} tests: {1} passed, {2} passed when run directly, {3} failed, {4} known pre-existing failures, {5} launch-blocked, {6} skipped, {7} timed out" -f `
        $tests.Count, (& $count 'Passed'), (& $count 'PassedDirect'), (& $count 'Failed'), (& $count 'KnownFailure'),
        (& $count 'LaunchBlocked'), (& $count 'Skipped'), (& $count 'Timeout')
    $notable = @($tests | Where-Object { $_.Category -notin @('Passed') } |
        ForEach-Object { "$($_.Name): $($_.Category) (ctest: $($_.Status)$(if ($_.Direct) { "; direct: $($_.Direct)" }))" })
    $ok = ($tests.Count -gt 0) -and ((& $count 'Failed') -eq 0) -and ((& $count 'Timeout') -eq 0)
    Add-Stage $Name $summary ($notable -join "`n") $ok
    return $tests
}

# --- Builds ------------------------------------------------------------------
$cudaCache = @(
    '-DOPTIMSOLVER_ENABLE_CUDA=ON', '-DOPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE=ON',
    '-DPDLP_BUILD_TESTS=ON', '-DQP_BUILD_TESTS=ON', '-DPDLP_BUILD_TOOLS=ON', '-DQP_BUILD_BENCH=ON',
    "-DCMAKE_CUDA_ARCHITECTURES=$CudaArchitectures")

if (-not $SkipBuild -and -not $SkipCpuBuild) {
    $logDir = New-OutputDirectory (Join-Path $OutDir 'build-cpu')
    $ok = Invoke-FreshBuild -SourceDir $root -BuildDir $CpuBuildDir -LogDir $logDir `
        -CacheArguments @('-DOPTIMSOLVER_ENABLE_CUDA=OFF', '-DPDLP_BUILD_TESTS=ON', '-DQP_BUILD_TESTS=ON')
    $warnings = Get-BuildWarnings (Join-Path $logDir 'build.log')
    $cache = Join-Path $CpuBuildDir 'CMakeCache.txt'
    $cudaInCache = (Test-Path $cache) -and (Select-String -Path $cache -Pattern '^CMAKE_CUDA_COMPILER:' -Quiet)
    Add-Stage 'CPU-only Release build' $(if ($ok) { "built; $($warnings.Count) distinct warning lines; CUDA compiler configured: $cudaInCache" } else { 'FAILED (see build-cpu logs)' }) `
        (($warnings | Select-Object -First 20) -join "`n") ($ok -and -not $cudaInCache)
}
if (-not $SkipBuild) {
    $logDir = New-OutputDirectory (Join-Path $OutDir 'build-cuda')
    $ok = Invoke-FreshBuild -SourceDir $root -BuildDir $CudaBuildDir -LogDir $logDir -CacheArguments $cudaCache
    $warnings = Get-BuildWarnings (Join-Path $logDir 'build.log')
    $cudaWarnings = @($warnings | Where-Object { $_ -match '\.cu[h]?\b' })
    $missing = @('pdlp_cuda_tests', 'qp_cuda_tests', 'optimsolver', 'pdlp_bench', 'qp_bench' |
        Where-Object { -not (Find-BuiltExecutable -BuildDir $CudaBuildDir -Name $_) })
    $result = if ($ok) { "built; $($warnings.Count) distinct warning lines ($($cudaWarnings.Count) from CUDA sources); missing targets: $(if ($missing) { $missing -join ', ' } else { 'none' })" } else { 'FAILED (see build-cuda logs)' }
    Add-Stage 'CUDA Release build' $result (($cudaWarnings + ($warnings | Select-Object -First 20)) -join "`n") ($ok -and $missing.Count -eq 0)
}

# --- Tests -------------------------------------------------------------------
if (-not $SkipCpuTests -and (Test-Path (Join-Path $CpuBuildDir 'CTestTestfile.cmake'))) {
    [void](Invoke-CTestStage -Name 'CTest (CPU-only build)' -BuildDir $CpuBuildDir)
}
[void](Invoke-CTestStage -Name 'CTest (CUDA build, full suite)' -BuildDir $CudaBuildDir)
$cudaTests = Invoke-CTestStage -Name 'CTest (CUDA build, -L cuda)' -BuildDir $CudaBuildDir -ExtraArgs @('-L', 'cuda')
if (@($cudaTests | Where-Object { $_.Category -eq 'Passed' }).Count -ne 2) {
    Add-Stage 'CUDA label coverage' 'expected exactly 2 passing CUDA-labelled tests' '' $false
}

foreach ($name in 'pdlp_cuda_tests', 'qp_cuda_tests') {
    $exe = Find-BuiltExecutable -BuildDir $CudaBuildDir -Name $name
    if (-not $exe) { Add-Stage "direct $name" 'not built' '' $false; continue }
    $run = Invoke-Native -FilePath $exe -WorkingDirectory (Split-Path $exe) -Environment @{ OPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE = '1' } `
        -LogPath (Join-Path $OutDir "$name.direct.log")
    $verdict = Get-RunVerdict $run
    $cases = ([regex]::Matches($run.Output, '(?m)^\s+pass\s')).Count
    $fails = ([regex]::Matches($run.Output, '(?m)^\s+FAIL\s')).Count
    Add-Stage "direct $name" "$verdict; $cases cases passed, $fails failed$(if ($run.LaunchError) { "; $($run.LaunchError)" })" '' ($verdict -eq 'Passed')
}

# --- Sanitizers and stability ------------------------------------------------
$powershell = (Get-Process -Id $PID).Path
if (-not $SkipSanitizers) {
    $sanitizerArgs = @('-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', (Join-Path $PSScriptRoot 'run_cuda_sanitizers.ps1'),
        '-BuildDir', $CudaBuildDir, '-OutDir', (Join-Path $OutDir 'sanitizers'))
    if ($Racecheck) { $sanitizerArgs += '-Racecheck' }
    $run = Invoke-Native -FilePath $powershell -Arguments $sanitizerArgs -LogPath (Join-Path $OutDir 'sanitizers.log')
    $lines = @($run.Output -split "`r?`n" | Where-Object { $_ -match '^(Passed|Failed|LaunchBlocked|MISSING)' })
    Add-Stage 'Compute Sanitizer' "exit $($run.ExitCode)" ($lines -join "`n") ($run.Launched -and $run.ExitCode -eq 0)
}
if (-not $SkipDeterminism) {
    $run = Invoke-Native -FilePath $powershell -LogPath (Join-Path $OutDir 'determinism.log') -Arguments @(
        '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', (Join-Path $PSScriptRoot 'check_cuda_determinism.ps1'),
        '-BuildDir', $CudaBuildDir, '-Runs', "$DeterminismRuns", '-OutDir', (Join-Path $OutDir 'determinism'))
    $lines = @($run.Output -split "`r?`n" | Where-Object { $_ -match 'runs passed' })
    Add-Stage "Determinism ($DeterminismRuns runs each)" "exit $($run.ExitCode)" ($lines -join "`n") ($run.Launched -and $run.ExitCode -eq 0)
}

# --- Summary -----------------------------------------------------------------
$hostInfo = Get-HostDescription
$commit = (Invoke-Native -FilePath (Get-Command git).Source -Arguments @('describe', '--always', '--dirty')).Output.Trim()
$md = @("# CUDA validation summary", '',
    "* Commit: $commit", "* CPU: $($hostInfo.Cpu) ($($hostInfo.LogicalProcessors) logical processors)",
    "* GPU: $($hostInfo.Gpu)", "* CUDA Toolkit: $($hostInfo.CudaToolkit)", "* OS: $($hostInfo.Os)",
    "* CUDA architectures: $CudaArchitectures", '', '| Stage | Result |', '| :--- | :--- |')
foreach ($key in $stages.Keys) { $md += "| $key | $(if ($stages[$key].Ok) { '' } else { '**FAIL** ' })$($stages[$key].Result) |" }
foreach ($key in $stages.Keys) {
    if ($stages[$key].Detail) { $md += @('', "## $key", '', '```', $stages[$key].Detail, '```') }
}
$md | Set-Content -Path (Join-Path $OutDir 'validation-summary.md') -Encoding UTF8
[pscustomobject]@{ Commit = $commit; Host = $hostInfo; Stages = $stages; Ok = $overallOk } |
    ConvertTo-Json -Depth 5 | Set-Content -Path (Join-Path $OutDir 'validation-summary.json') -Encoding UTF8
Write-Host ''
Write-Host "Summary: $(Join-Path $OutDir 'validation-summary.md')"
Write-Host $(if ($overallOk) { 'VALIDATION PASSED' } else { 'VALIDATION FAILED' })
if ($overallOk) { exit 0 } else { exit 1 }

} finally {
    Pop-Location
}

<#
.SYNOPSIS
    Runs NVIDIA Compute Sanitizer over the CUDA test binaries.

.DESCRIPTION
    Runs memcheck and initcheck (and, with -Racecheck, racecheck) on
    pdlp_cuda_tests and qp_cuda_tests from an existing CUDA build. A check
    passes only if the test binary itself passes and the sanitizer reports
    "ERROR SUMMARY: 0 errors". Full logs go to -OutDir; a summary is printed
    and written as sanitizers.json.

    Exit code: 0 all passed, 1 a sanitizer or test failure, 2 a setup problem
    (no sanitizer, missing binary, launch blocked by the OS).

.EXAMPLE
    tools/cuda/run_cuda_sanitizers.ps1 -BuildDir build-cuda
#>
param(
    [Parameter(Mandatory)] [string] $BuildDir,
    [string] $OutDir = (Join-Path $BuildDir 'validation/sanitizers'),
    [string[]] $Tools = @('memcheck', 'initcheck'),
    [switch] $Racecheck,
    [string[]] $Binaries = @('pdlp_cuda_tests', 'qp_cuda_tests')
)

$ErrorActionPreference = 'Stop'
Import-Module (Join-Path $PSScriptRoot 'CudaValidation.psm1') -Force

if ($Racecheck -and $Tools -notcontains 'racecheck') { $Tools += 'racecheck' }
$OutDir = New-OutputDirectory $OutDir
$sanitizer = Find-ComputeSanitizer
if (-not $sanitizer) {
    Write-Error 'compute-sanitizer not found (CUDA Toolkit bin directory on PATH, or CUDA_PATH set?)' -ErrorAction Continue
    exit 2
}
Write-Host "compute-sanitizer: $sanitizer"

$results = @()
$setupProblem = $false
foreach ($name in $Binaries) {
    $exe = Find-BuiltExecutable -BuildDir $BuildDir -Name $name
    if (-not $exe) {
        Write-Host "MISSING  $name (not built under $BuildDir)"
        $results += [pscustomobject]@{ Binary = $name; Tool = '-'; Verdict = 'Missing'; Errors = $null; TestExit = $null; Seconds = 0 }
        $setupProblem = $true
        continue
    }
    foreach ($tool in $Tools) {
        $log = Join-Path $OutDir "$name.$tool.log"
        $sanitizerArgs = @('--tool', $tool, '--error-exitcode', '99')
        if ($tool -eq 'memcheck') { $sanitizerArgs += @('--leak-check', 'full') }
        $run = Invoke-Native -FilePath $sanitizer -Arguments ($sanitizerArgs + @($exe)) `
            -WorkingDirectory (Split-Path $exe) -LogPath $log `
            -Environment @{ OPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE = '1' }
        $errors = $null
        if ($run.Output -match 'ERROR SUMMARY:\s+(\d+) error') { $errors = [int]$Matches[1] }
        $testPassed = $run.Output -match 'All (PDLP|QP) CUDA tests passed'
        $verdict = if (-not $run.Launched) { 'LaunchBlocked' }
                   elseif ($errors -eq 0 -and $testPassed -and $run.ExitCode -eq 0) { 'Passed' }
                   else { 'Failed' }
        if ($verdict -eq 'LaunchBlocked') { $setupProblem = $true }
        $results += [pscustomobject]@{
            Binary = $name; Tool = $tool; Verdict = $verdict; Errors = $errors
            TestExit = $run.ExitCode; TestsPassed = $testPassed; Seconds = [math]::Round($run.Seconds, 1)
            LaunchError = $run.LaunchError
        }
        $summary = if ($errors -ne $null) { "ERROR SUMMARY: $errors errors" } else { 'no ERROR SUMMARY line' }
        Write-Host ("{0,-8} {1,-16} {2,-10} {3}  ({4} s)" -f $verdict, $name, $tool, $summary, [math]::Round($run.Seconds, 1))
        if ($run.LaunchError) { Write-Host "         $($run.LaunchError)" }
    }
}

$results | ConvertTo-Json -Depth 3 | Set-Content -Path (Join-Path $OutDir 'sanitizers.json') -Encoding UTF8
if (@($results | Where-Object { $_.Verdict -eq 'Failed' }).Count -gt 0) { exit 1 }
if ($setupProblem) { exit 2 }
exit 0

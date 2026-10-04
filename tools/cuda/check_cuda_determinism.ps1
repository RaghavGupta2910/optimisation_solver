<#
.SYNOPSIS
    Repeats the CUDA test binaries and checks that every run passes and
    reports the same deterministic results.

.DESCRIPTION
    Runs pdlp_cuda_tests and qp_cuda_tests -Runs times each (default 10). Each
    run must exit 0. From every run's output the script keeps the per-test
    pass/fail lines and the reported CPU and CUDA iteration counts and transfer
    byte counts, drops wall-clock figures (setup and KKT seconds), and requires
    the result to be identical across all runs. Timings are never compared.

    The binaries themselves also assert bitwise run-to-run reproducibility of
    CUDA solves ("cuda solve is reproducible" cases), so this script checks
    reproducibility both within a process and across processes.

    Exit code: 0 stable, 1 a failure or a difference, 2 a setup problem.

.EXAMPLE
    tools/cuda/check_cuda_determinism.ps1 -BuildDir build-cuda -Runs 10
#>
param(
    [Parameter(Mandatory)] [string] $BuildDir,
    [int] $Runs = 10,
    [string] $OutDir = (Join-Path $BuildDir 'validation/determinism'),
    [string[]] $Binaries = @('pdlp_cuda_tests', 'qp_cuda_tests')
)

$ErrorActionPreference = 'Stop'
Import-Module (Join-Path $PSScriptRoot 'CudaValidation.psm1') -Force
$OutDir = New-OutputDirectory $OutDir

# Keeps only what must not change between runs.
function Get-Fingerprint([string] $Text) {
    $lines = foreach ($line in ($Text -split "`r?`n")) {
        if ($line -match '^\s+(pass|FAIL)\s+') {
            $line.Trim()
        } elseif ($line -match 'iterations cpu') {
            $l = $line -replace '(setup|kkt)\s+[-+0-9.eE]+\s*s', '$1 <t>'
            $l.Trim()
        }
    }
    return ($lines -join "`n")
}

$results = @()
$problem = $false
$setupProblem = $false
foreach ($name in $Binaries) {
    $exe = Find-BuiltExecutable -BuildDir $BuildDir -Name $name
    if (-not $exe) {
        Write-Host "MISSING  $name"
        $setupProblem = $true
        continue
    }
    $reference = $null
    $passed = 0; $blocked = 0; $failed = 0; $differing = 0
    for ($i = 1; $i -le $Runs; $i++) {
        $run = Invoke-Native -FilePath $exe -WorkingDirectory (Split-Path $exe) `
            -LogPath (Join-Path $OutDir "$name.run$i.log") `
            -Environment @{ OPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE = '1' }
        $verdict = Get-RunVerdict $run
        switch ($verdict) {
            'Passed' { $passed++ }
            'LaunchBlocked' { $blocked++ }
            default { $failed++ }
        }
        if ($verdict -ne 'Passed') {
            Write-Host "  $name run ${i}: $verdict $($run.LaunchError)"
            continue
        }
        $fingerprint = Get-Fingerprint $run.Output
        if ($null -eq $reference) {
            $reference = $fingerprint
            Set-Content -Path (Join-Path $OutDir "$name.fingerprint.txt") -Value $fingerprint -Encoding UTF8
        } elseif ($fingerprint -ne $reference) {
            $differing++
            Set-Content -Path (Join-Path $OutDir "$name.run$i.fingerprint.txt") -Value $fingerprint -Encoding UTF8
            Write-Host "  $name run ${i}: deterministic output differs from run 1"
        }
    }
    $checks = if ($reference) { @($reference -split "`n" | Where-Object { $_ -match '^pass' }).Count } else { 0 }
    $ok = ($passed -eq $Runs -and $differing -eq 0)
    if ($failed -gt 0 -or $differing -gt 0) { $problem = $true }
    if ($blocked -gt 0) { $setupProblem = $true }
    $results += [pscustomobject]@{
        Binary = $name; Runs = $Runs; Passed = $passed; Failed = $failed; LaunchBlocked = $blocked
        Differing = $differing; CasesPerRun = $checks; Stable = $ok
    }
    Write-Host ("{0,-16} {1}/{2} runs passed, {3} cases per run, {4} differing fingerprints -> {5}" -f
        $name, $passed, $Runs, $checks, $differing, $(if ($ok) { 'STABLE' } else { 'NOT STABLE' }))
}

$results | ConvertTo-Json -Depth 3 | Set-Content -Path (Join-Path $OutDir 'determinism.json') -Encoding UTF8
if ($problem) { exit 1 }
if ($setupProblem) { exit 2 }
exit 0

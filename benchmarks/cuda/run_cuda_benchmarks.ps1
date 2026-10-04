<#
.SYNOPSIS
    CPU vs CUDA end-to-end benchmarks for PDLP and for the hybrid QP backend.

.DESCRIPTION
    Drives pdlp_bench and `qp_bench sweep` from a Release CUDA build (configured
    with -DPDLP_BUILD_TOOLS=ON -DQP_BUILD_BENCH=ON) over a fixed list of
    generated instances. Every instance is generated from a fixed seed, so the
    same command reproduces the same problems.

    For every instance and backend the harness performs one untimed warm-up
    solve and -Runs timed solves (default 5) of the complete solver call, and
    reports the median wall time. For CUDA that time includes device setup,
    host-to-device upload, all kernels and reductions, host-side checks and
    the final download. CPU and CUDA are run back to back per instance.

    Outputs in -OutDir:
      pdlp.json / qp.json     every raw result record (all timed runs)
      pdlp.csv  / qp.csv      one row per instance, CPU next to CUDA
      summary.md              Markdown tables and the host description

    Speedup = CPU median / CUDA median. Objective error = |CPU obj - CUDA obj|.

.EXAMPLE
    benchmarks/cuda/run_cuda_benchmarks.ps1 -BuildDir build-validate-cuda
#>
param(
    [Parameter(Mandatory)] [string] $BuildDir,
    [string] $OutDir = (Join-Path 'benchmarks/results/cuda' (Get-Date -Format 'yyyy-MM-dd')),
    [int] $Runs = 5,
    [ValidateSet('pdlp', 'qp')] [string[]] $Engines = @('pdlp', 'qp'),
    # rows x columns x nonzeros-per-row; columns = 2 * rows, so nnz ~ 10 * rows.
    [string[]] $PdlpInstances = @(
        '1000x2000x10', '10000x20000x10', '25000x50000x10', '50000x100000x10', '75000x150000x10',
        '100000x200000x10', '150000x300000x10', '200000x400000x10', '500000x1000000x10'),
    [int] $PdlpIterationLimit = 200000,
    # variables (= constraint rows) x nonzeros per constraint row
    [string[]] $QpInstances = @('2000x8', '5000x8', '10000x8', '20000x8', '30000x8'),
    [int] $Threads = 0,
    [int] $TimeoutSeconds = 7200,
    # Reuse a completed result from an earlier, interrupted run into the same
    # -OutDir (same backend, thread count and number of runs) instead of
    # re-measuring it.
    [switch] $Resume
)

$ErrorActionPreference = 'Stop'
# Tables and CSV must not depend on the machine's locale (digit grouping,
# decimal separator).
[System.Threading.Thread]::CurrentThread.CurrentCulture = [System.Globalization.CultureInfo]::InvariantCulture
Import-Module (Join-Path $PSScriptRoot '..\..\tools\cuda\CudaValidation.psm1') -Force
$OutDir = New-OutputDirectory $OutDir
$logDir = New-OutputDirectory (Join-Path $OutDir 'logs')

function Read-BenchLog([string] $Log) {
    if (-not (Test-Path $Log)) { return $null }
    $line = (Get-Content $Log | Where-Object { $_ -like 'json *' } | Select-Object -Last 1)
    if (-not $line) { return $null }
    return ($line.Substring(5) | ConvertFrom-Json)
}

function Invoke-Bench([string] $Exe, [string[]] $Arguments, [string] $Log, [string] $Backend) {
    if ($Resume) {
        $previous = Read-BenchLog $Log
        if ($previous -and $previous.requested_backend -eq $Backend -and $previous.threads -eq $Threads -and
            @($previous.run_seconds).Count -eq $Runs) {
            Write-Host "  (reusing $Log)"
            return $previous
        }
    }
    $run = Invoke-Native -FilePath $Exe -Arguments $Arguments -WorkingDirectory (Split-Path $Exe) `
        -LogPath $Log -TimeoutSeconds $TimeoutSeconds
    $verdict = Get-RunVerdict $run
    if ($verdict -ne 'Passed') {
        throw "$(Split-Path -Leaf $Exe) $($Arguments -join ' '): $verdict $($run.LaunchError) (log: $Log)"
    }
    $line = ($run.Output -split "`r?`n" | Where-Object { $_ -like 'json *' } | Select-Object -Last 1)
    if (-not $line) { throw "no json result line in $Log" }
    return ($line.Substring(5) | ConvertFrom-Json)
}

function Get-Spread($Record) {
    $sorted = @($Record.run_seconds | Sort-Object)
    return [pscustomobject]@{ Min = $sorted[0]; Max = $sorted[-1] }
}

function Format-Seconds($Value) {
    if ($null -eq $Value) { return 'n/a' }
    if ($Value -lt 0.1) { return ('{0:N4}' -f $Value) }
    return ('{0:N3}' -f $Value)
}

function Get-Comparison($Cpu, $Gpu, [string] $Name) {
    $objectiveError = $null
    $relativeError = $null
    if ($null -ne $Cpu.objective -and $null -ne $Gpu.objective) {
        $objectiveError = [math]::Abs($Cpu.objective - $Gpu.objective)
        $relativeError = $objectiveError / [math]::Max(1.0, [math]::Abs($Cpu.objective))
    }
    $cpuSpread = Get-Spread $Cpu
    $gpuSpread = Get-Spread $Gpu
    return [pscustomobject][ordered]@{
        instance = $Name; rows = $Cpu.rows; columns = $Cpu.columns; nonzeros = $Cpu.nonzeros
        cpu_status = $Cpu.status; cuda_status = $Gpu.status
        cpu_iterations = $Cpu.iterations; cuda_iterations = $Gpu.iterations
        cuda_executed = $Gpu.executed_backend
        cpu_seconds = $Cpu.median_seconds; cuda_seconds = $Gpu.median_seconds
        cpu_min_seconds = $cpuSpread.Min; cpu_max_seconds = $cpuSpread.Max
        cuda_min_seconds = $gpuSpread.Min; cuda_max_seconds = $gpuSpread.Max
        speedup = $Cpu.median_seconds / $Gpu.median_seconds
        cuda_setup_seconds = $Gpu.setup_seconds
        cuda_h2d_bytes = $Gpu.h2d_bytes; cuda_d2h_bytes = $Gpu.d2h_bytes
        cpu_kkt_solve_seconds = $(if ($Cpu.PSObject.Properties['kkt_solve_seconds']) { $Cpu.kkt_solve_seconds } else { $null })
        cuda_kkt_solve_seconds = $(if ($Gpu.PSObject.Properties['kkt_solve_seconds']) { $Gpu.kkt_solve_seconds } else { $null })
        cuda_host_check_seconds = $(if ($Gpu.PSObject.Properties['host_check_seconds']) { $Gpu.host_check_seconds } else { $null })
        cpu_objective = $Cpu.objective; cuda_objective = $Gpu.objective
        objective_error = $objectiveError; relative_objective_error = $relativeError
    }
}

$hostInfo = Get-HostDescription
$commit = (Invoke-Native -FilePath (Get-Command git).Source -Arguments @('describe', '--always', '--dirty')).Output.Trim()
Write-Host "CPU: $($hostInfo.Cpu); GPU: $($hostInfo.Gpu); CUDA $($hostInfo.CudaToolkit); commit $commit"
$md = @('# CPU vs CUDA benchmarks', '',
    "* Commit: $commit",
    "* CPU: $($hostInfo.Cpu) ($($hostInfo.LogicalProcessors) logical processors; solver threads: $(if ($Threads -eq 0) { 'all' } else { $Threads }))",
    "* GPU: $($hostInfo.Gpu)", "* CUDA Toolkit: $($hostInfo.CudaToolkit)", "* OS: $($hostInfo.Os)",
    "* Method: Release build, double precision, 1 untimed warm-up + $Runs timed end-to-end solves per backend; median reported; CPU and CUDA run back to back per instance",
    '* Speedup = CPU median / CUDA median; objective error = abs(CPU objective - CUDA objective)')

if ($Engines -contains 'pdlp') {
    $exe = Find-BuiltExecutable -BuildDir $BuildDir -Name 'pdlp_bench'
    if (-not $exe) { throw "pdlp_bench not found under $BuildDir (configure with -DPDLP_BUILD_TOOLS=ON)" }
    $raw = @(); $rows = @()
    foreach ($instance in $PdlpInstances) {
        $dims = $instance -split 'x'
        $records = @{}
        foreach ($backend in 'cpu', 'cuda') {
            Write-Host "pdlp $instance $backend ..."
            $records[$backend] = Invoke-Bench $exe @($dims[0], $dims[1], $dims[2], "$PdlpIterationLimit", "$Threads", $backend, "$Runs") `
                (Join-Path $logDir "pdlp_${instance}_$backend.log") $backend
            $raw += $records[$backend]
        }
        $row = Get-Comparison $records['cpu'] $records['cuda'] $instance
        $rows += $row
        Write-Host ("  nnz {0}: cpu {1} s ({2}, {3} it), cuda {4} s ({5}, {6} it), speedup {7:N2}x, obj err {8}" -f
            $row.nonzeros, (Format-Seconds $row.cpu_seconds), $row.cpu_status, $row.cpu_iterations,
            (Format-Seconds $row.cuda_seconds), $row.cuda_status, $row.cuda_iterations, $row.speedup, $row.objective_error)
    }
    $raw | ConvertTo-Json -Depth 4 | Set-Content (Join-Path $OutDir 'pdlp.json') -Encoding UTF8
    $rows | Export-Csv -NoTypeInformation -Path (Join-Path $OutDir 'pdlp.csv') -Encoding UTF8
    $md += @('', '## PDLP: CPU vs CUDA', '',
        "Generated feasible LPs (pdlp_bench, seed 12345): rows x 2*rows columns, 10 nonzeros per row, iteration limit $PdlpIterationLimit, default tolerances (1e-6), polishing off.", '',
        '| Instance | NNZ | CPU (s) | CUDA (s) | Speedup | Objective Error | CPU / CUDA iterations | CPU / CUDA status | CUDA setup (s) | CPU min-max (s) | CUDA min-max (s) |',
        '| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: |')
    foreach ($r in $rows) {
        $md += ('| {0} | {1:N0} | {2} | {3} | {4:N2}x | {5:0.00e+00} | {6} / {7} | {8} / {9} | {10} | {11}-{12} | {13}-{14} |' -f
            $r.instance, $r.nonzeros, (Format-Seconds $r.cpu_seconds), (Format-Seconds $r.cuda_seconds), $r.speedup,
            $r.objective_error, $r.cpu_iterations, $r.cuda_iterations, $r.cpu_status, $r.cuda_status,
            (Format-Seconds $r.cuda_setup_seconds), (Format-Seconds $r.cpu_min_seconds), (Format-Seconds $r.cpu_max_seconds),
            (Format-Seconds $r.cuda_min_seconds), (Format-Seconds $r.cuda_max_seconds))
    }
}

if ($Engines -contains 'qp') {
    $exe = Find-BuiltExecutable -BuildDir $BuildDir -Name 'qp_bench'
    if (-not $exe) { throw "qp_bench not found under $BuildDir (configure with -DQP_BUILD_BENCH=ON)" }
    $raw = @(); $rows = @()
    foreach ($instance in $QpInstances) {
        $dims = $instance -split 'x'
        $records = @{}
        foreach ($backend in 'cpu', 'cuda') {
            Write-Host "qp $instance $backend ..."
            $records[$backend] = Invoke-Bench $exe @('sweep', $backend, $dims[0], $dims[1], "$Runs", "$Threads") `
                (Join-Path $logDir "qp_${instance}_$backend.log") $backend
            $raw += $records[$backend]
        }
        $row = Get-Comparison $records['cpu'] $records['cuda'] $instance
        $rows += $row
        Write-Host ("  nnz {0}: cpu {1} s ({2}), hybrid cuda {3} s ({4}), speedup {5:N2}x, cuda kkt {6} s, obj err {7}" -f
            $row.nonzeros, (Format-Seconds $row.cpu_seconds), $row.cpu_status, (Format-Seconds $row.cuda_seconds),
            $row.cuda_status, $row.speedup, (Format-Seconds $row.cuda_kkt_solve_seconds), $row.objective_error)
    }
    $raw | ConvertTo-Json -Depth 4 | Set-Content (Join-Path $OutDir 'qp.json') -Encoding UTF8
    $rows | Export-Csv -NoTypeInformation -Path (Join-Path $OutDir 'qp.csv') -Encoding UTF8
    $md += @('', '## QP: CPU vs hybrid CUDA', '',
        'Generated banded convex QPs (`qp_bench sweep`, seed 2024): n variables, n constraint rows with 8 nonzeros each, tridiagonal P, default ADMM options. The CUDA backend is hybrid: the KKT factorisation and solves stay on the CPU.', '',
        '| Instance | NNZ (P+A) | CPU (s) | Hybrid CUDA (s) | Speedup | CUDA KKT solve (s) | CUDA setup (s) | H2D / D2H (MB) | CPU / CUDA iterations | CPU / CUDA status | Objective Error |',
        '| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: |')
    foreach ($r in $rows) {
        $md += ('| {0} | {1:N0} | {2} | {3} | {4:N2}x | {5} | {6} | {7:N1} / {8:N1} | {9} / {10} | {11} / {12} | {13:0.00e+00} |' -f
            $r.instance, $r.nonzeros, (Format-Seconds $r.cpu_seconds), (Format-Seconds $r.cuda_seconds), $r.speedup,
            (Format-Seconds $r.cuda_kkt_solve_seconds), (Format-Seconds $r.cuda_setup_seconds),
            ($r.cuda_h2d_bytes / 1MB), ($r.cuda_d2h_bytes / 1MB), $r.cpu_iterations, $r.cuda_iterations,
            $r.cpu_status, $r.cuda_status, $r.objective_error)
    }
}

$md | Set-Content (Join-Path $OutDir 'summary.md') -Encoding UTF8
Write-Host "Results: $OutDir"

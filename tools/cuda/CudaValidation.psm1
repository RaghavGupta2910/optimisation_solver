# Shared helpers for the CUDA validation and benchmark scripts.
#
# Works on Windows PowerShell 5.1 and PowerShell 7 (Windows or Linux).
#
# Every child process is started through System.Diagnostics.Process rather than
# the call operator, so the result is unambiguous: either the process launched
# and has an exit code, or it never launched (for example because Windows
# Application Control blocked a freshly built executable) and there is no exit
# code at all. $LASTEXITCODE is never consulted, so a stale value from an
# earlier command cannot be mistaken for a result.

Set-StrictMode -Version 3.0

$script:IsWindowsHost = ($env:OS -eq 'Windows_NT')

function Get-RepoRoot {
    return (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
}

# Runs an executable and returns
#   Launched   $true if the process started
#   ExitCode   its exit code ($null if it did not start)
#   Output     combined stdout + stderr text
#   LaunchError the reason it could not start ($null if it did)
#   Seconds    wall time
# Output is also written to LogPath when given. Never throws for a non-zero
# exit or a launch failure; callers decide what each means.
function Invoke-Native {
    param(
        [Parameter(Mandatory)] [string] $FilePath,
        [string[]] $Arguments = @(),
        [string] $WorkingDirectory = (Get-Location).Path,
        [hashtable] $Environment = @{},
        [string] $LogPath,
        [int] $TimeoutSeconds = 0
    )
    $info = New-Object System.Diagnostics.ProcessStartInfo
    $info.FileName = $FilePath
    $info.Arguments = (($Arguments | ForEach-Object { ConvertTo-ProcessArgument $_ }) -join ' ')
    $info.WorkingDirectory = $WorkingDirectory
    $info.UseShellExecute = $false
    $info.RedirectStandardOutput = $true
    $info.RedirectStandardError = $true
    $info.CreateNoWindow = $true
    foreach ($key in $Environment.Keys) {
        $info.EnvironmentVariables[$key] = [string]$Environment[$key]
    }

    $process = New-Object System.Diagnostics.Process
    $process.StartInfo = $info
    $watch = [System.Diagnostics.Stopwatch]::StartNew()
    try {
        [void]$process.Start()
    } catch {
        $inner = $_.Exception
        while ($inner.InnerException) { $inner = $inner.InnerException }
        $message = $inner.Message
        if ($inner -is [System.ComponentModel.Win32Exception]) {
            $message = "$message (Win32 error $($inner.NativeErrorCode))"
        }
        if ($LogPath) { Set-Content -Path $LogPath -Value "LAUNCH FAILED: $message" -Encoding UTF8 }
        return [pscustomobject]@{
            Launched = $false; ExitCode = $null; Output = ''; LaunchError = $message
            Seconds = $watch.Elapsed.TotalSeconds; TimedOut = $false
        }
    }
    $stdout = $process.StandardOutput.ReadToEndAsync()
    $stderr = $process.StandardError.ReadToEndAsync()
    $timedOut = $false
    if ($TimeoutSeconds -gt 0) {
        if (-not $process.WaitForExit($TimeoutSeconds * 1000)) {
            $timedOut = $true
            try { $process.Kill() } catch { }
        }
    }
    $process.WaitForExit()
    $watch.Stop()
    $text = $stdout.Result
    if ($stderr.Result) { $text = $text + $stderr.Result }
    if ($LogPath) { Set-Content -Path $LogPath -Value $text -Encoding UTF8 }
    return [pscustomobject]@{
        Launched = $true; ExitCode = $process.ExitCode; Output = $text; LaunchError = $null
        Seconds = $watch.Elapsed.TotalSeconds; TimedOut = $timedOut
    }
}

function ConvertTo-ProcessArgument([string] $Value) {
    if ($Value -eq '') { return '""' }
    if ($Value -notmatch '[\s"]') { return $Value }
    # Windows command-line quoting: double backslashes that precede a quote.
    $escaped = [regex]::Replace($Value, '(\\*)"', '$1$1\"')
    $escaped = [regex]::Replace($escaped, '(\\+)$', '$1$1')
    return '"' + $escaped + '"'
}

# Classifies an Invoke-Native result for reporting.
function Get-RunVerdict($Run, [int[]] $SkipCodes = @()) {
    if (-not $Run.Launched) { return 'LaunchBlocked' }
    if ($Run.TimedOut) { return 'Timeout' }
    if ($Run.ExitCode -eq 0) { return 'Passed' }
    if ($SkipCodes -contains $Run.ExitCode) { return 'Skipped' }
    return 'Failed'
}

# Finds a built executable by name anywhere under the build tree, preferring a
# Release configuration directory for multi-config generators.
function Find-BuiltExecutable {
    param([Parameter(Mandatory)] [string] $BuildDir, [Parameter(Mandatory)] [string] $Name)
    $file = if ($script:IsWindowsHost) { "$Name.exe" } else { $Name }
    $found = @(Get-ChildItem -Path $BuildDir -Recurse -File -Filter $file -ErrorAction SilentlyContinue |
        Where-Object { $_.FullName -notmatch '[\\/]CMakeFiles[\\/]' })
    if ($found.Count -eq 0) { return $null }
    $release = @($found | Where-Object { $_.FullName -match '[\\/]Release[\\/]' })
    if ($release.Count -gt 0) { return $release[0].FullName }
    return $found[0].FullName
}

function Find-ComputeSanitizer {
    $candidates = @()
    $command = Get-Command compute-sanitizer -ErrorAction SilentlyContinue
    if ($command) {
        $dir = Split-Path $command.Source
        $candidates += (Join-Path $dir '..\compute-sanitizer\compute-sanitizer.exe')
        $candidates += (Join-Path $dir '..\compute-sanitizer\compute-sanitizer')
        $candidates += $command.Source
    }
    foreach ($root in @($env:CUDA_PATH, '/usr/local/cuda')) {
        if ($root) {
            $candidates += (Join-Path $root 'compute-sanitizer\compute-sanitizer.exe')
            $candidates += (Join-Path $root 'compute-sanitizer/compute-sanitizer')
            $candidates += (Join-Path $root 'bin/compute-sanitizer')
        }
    }
    foreach ($candidate in $candidates) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    return $null
}

# Configures and builds a fresh tree. Returns $true on success.
function Invoke-FreshBuild {
    param(
        [Parameter(Mandatory)] [string] $SourceDir,
        [Parameter(Mandatory)] [string] $BuildDir,
        [string[]] $CacheArguments = @(),
        [string] $Config = 'Release',
        [Parameter(Mandatory)] [string] $LogDir,
        [switch] $KeepExisting
    )
    if ((Test-Path $BuildDir) -and -not $KeepExisting) {
        Remove-Item -Recurse -Force $BuildDir
    }
    $cmake = (Get-Command cmake -ErrorAction Stop).Source
    $configure = Invoke-Native -FilePath $cmake -LogPath (Join-Path $LogDir 'configure.log') -Arguments (
        @('-S', $SourceDir, '-B', $BuildDir, "-DCMAKE_BUILD_TYPE=$Config") + $CacheArguments)
    if ((Get-RunVerdict $configure) -ne 'Passed') { return $false }
    $parallel = [Environment]::ProcessorCount
    $build = Invoke-Native -FilePath $cmake -LogPath (Join-Path $LogDir 'build.log') -Arguments @(
        '--build', $BuildDir, '--config', $Config, '--parallel', "$parallel")
    return ((Get-RunVerdict $build) -eq 'Passed')
}

# Warning lines from a build log, de-duplicated (MSBuild repeats each one).
function Get-BuildWarnings([string] $LogPath) {
    if (-not (Test-Path $LogPath)) { return @() }
    return @(Get-Content $LogPath | Where-Object { $_ -match '(?i)\bwarning\b' -and $_ -notmatch '0 Warning\(s\)' } |
        ForEach-Object { ($_ -replace '\s+\[[^\]]+\]$', '').Trim() } | Sort-Object -Unique)
}

function Get-HostDescription {
    $cpu = $null
    $gpu = $null
    try {
        if ($script:IsWindowsHost) {
            $cpu = (Get-CimInstance Win32_Processor | Select-Object -First 1).Name.Trim()
        } elseif (Test-Path /proc/cpuinfo) {
            $cpu = ((Get-Content /proc/cpuinfo | Where-Object { $_ -match '^model name' } | Select-Object -First 1) -split ':', 2)[1].Trim()
        }
    } catch { }
    $smi = Get-Command nvidia-smi -ErrorAction SilentlyContinue
    if ($smi) {
        $q = Invoke-Native -FilePath $smi.Source -Arguments @('--query-gpu=name,compute_cap,driver_version,memory.total', '--format=csv,noheader')
        if ($q.Launched -and $q.ExitCode -eq 0) { $gpu = ($q.Output -split "`n")[0].Trim() }
    }
    $nvcc = Get-Command nvcc -ErrorAction SilentlyContinue
    $toolkit = $null
    if ($nvcc) {
        $v = Invoke-Native -FilePath $nvcc.Source -Arguments @('--version')
        if ($v.Output -match 'release ([0-9.]+)') { $toolkit = $Matches[1] }
    }
    return [pscustomobject]@{
        Cpu = $cpu; LogicalProcessors = [Environment]::ProcessorCount; Gpu = $gpu; CudaToolkit = $toolkit
        Os = [System.Environment]::OSVersion.VersionString
    }
}

function New-OutputDirectory([string] $Path) {
    New-Item -ItemType Directory -Force -Path $Path | Out-Null
    return (Resolve-Path $Path).Path
}

Export-ModuleMember -Function Get-RepoRoot, Invoke-Native, Get-RunVerdict, Find-BuiltExecutable,
    Find-ComputeSanitizer, Invoke-FreshBuild, Get-BuildWarnings, Get-HostDescription, New-OutputDirectory

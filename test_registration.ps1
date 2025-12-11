param(
    [Parameter(Position=0)] [string]$fixed,
    [Parameter(Position=1)] [string]$moving,
    [Parameter(Position=2)] [string]$output,
    [string]$config,
    [string]$initial,
    [string]$fixedMask,
    [string]$transform,
    [double]$samplingPercentage,
    [switch]$DryRun
)

# 测试配准程序的示例脚本
# 使用方法: .\test_registration.ps1

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Test MI Registration" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# 检查可执行文件是否存在
$exe = ".\build\bin\Release\MIRegistration.exe"
if (-not (Test-Path $exe)) {
    Write-Host "`n[!] Error: Executable not found!" -ForegroundColor Red
    Write-Host "    Please build the project first using: .\build.ps1" -ForegroundColor Yellow
    exit 1
}

if (-not $fixed -or -not $moving -or -not $output) {
    Write-Host "`nPlease provide the following image files:" -ForegroundColor Yellow
    Write-Host "(Supported formats: .nrrd, .nii, .nii.gz, .mha, .mhd, DICOM)" -ForegroundColor Gray
    if (-not $fixed) { $fixed = Read-Host "Fixed image path" }
    if (-not $moving) { $moving = Read-Host "Moving image path" }
    if (-not $output) { $output = Read-Host "Output folder path (for .h5 transform file)" }
    if (-not $config) { $config = Read-Host "Optional: config json path (press Enter to use default)" }
    if (-not $initial) { $initial = Read-Host "Optional: initial transform .h5 file (press Enter to skip)" }
    if (-not $fixedMask) { $fixedMask = Read-Host "Optional: fixed mask .nrrd file (press Enter to skip)" }
    if (-not $transform) { $transform = Read-Host "Transform type (Rigid/Affine) [Rigid]"; if (-not $transform) { $transform = "Rigid" } }
}

# 检查输入文件是否存在
if (-not (Test-Path $fixed)) {
    Write-Host "`n[!] Error: Fixed image not found: $fixed" -ForegroundColor Red
    exit 1
}

if (-not (Test-Path $moving)) {
    Write-Host "`n[!] Error: Moving image not found: $moving" -ForegroundColor Red
    exit 1
}

# 运行配准
Write-Host "`n[*] Starting registration..." -ForegroundColor Green
Write-Host "    Fixed:  $fixed" -ForegroundColor Gray
Write-Host "    Moving: $moving" -ForegroundColor Gray
Write-Host "    Output: $output" -ForegroundColor Gray
Write-Host ""

# Ensure output folder exists
if (-not (Test-Path $output)) {
    New-Item -ItemType Directory -Path $output -Force | Out-Null
}

# Build argument list for executable
$argsList = @()
# Pass positional arguments in order (main.cpp supports positional 1-3)
$argsList += $fixed
$argsList += $moving
$argsList += $output
if ($config) { $argsList += "--config"; $argsList += $config }
if ($initial) { $argsList += "--initial"; $argsList += $initial }
# Pass fixed mask for local registration (only ROI voxels used)
if ($fixedMask) { $argsList += "--fixed-mask"; $argsList += $fixedMask }
# Only pass transform if explicitly provided by user
if ($PSBoundParameters.ContainsKey('transform')) { $argsList += "--transform"; $argsList += $transform }
if ($PSBoundParameters.ContainsKey('Verbose')) { $argsList += "--verbose" }
# Only pass samplingPercentage if explicitly provided by user (to not override config file)
if ($PSBoundParameters.ContainsKey('samplingPercentage')) { $argsList += "--sampling-percentage"; $argsList += $samplingPercentage }

Write-Host "Invoking: $exe $($argsList -join ' ')`n" -ForegroundColor DarkCyan
if (-not $DryRun) {
    & $exe @argsList
} else {
    Write-Host "DryRun: Skipping execution (DryRun enabled)" -ForegroundColor Yellow
}

if ($LASTEXITCODE -eq 0) {
    Write-Host "`n========================================" -ForegroundColor Cyan
    Write-Host "  Registration Completed!" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host "`nTransform (.h5) saved to folder: $output" -ForegroundColor Yellow
} else {
    Write-Host "`n[!] Registration failed with error code: $LASTEXITCODE" -ForegroundColor Red
}

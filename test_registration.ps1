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

# 获取输入文件
Write-Host "`nPlease provide the following image files:" -ForegroundColor Yellow
Write-Host "(Supported formats: .nrrd, .nii, .nii.gz, .mha, .mhd, DICOM)" -ForegroundColor Gray
$fixed = Read-Host "Fixed image path"
$moving = Read-Host "Moving image path"
$output = Read-Host "Output folder path (for .h5 transform file)"

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

& $exe $fixed $moving $output

if ($LASTEXITCODE -eq 0) {
    Write-Host "`n========================================" -ForegroundColor Cyan
    Write-Host "  Registration Completed!" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Cyan
    Write-Host "`nTransform (.h5) saved to folder: $output" -ForegroundColor Yellow
} else {
    Write-Host "`n[!] Registration failed with error code: $LASTEXITCODE" -ForegroundColor Red
}

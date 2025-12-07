# 构建脚本 - Windows PowerShell (使用vcpkg)
# 使用方法: .\build.ps1

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  MI Registration Build Script (vcpkg)" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

# 检查vcpkg环境变量
if (-not $env:VCPKG_ROOT) {
    Write-Host "`n[!] Warning: VCPKG_ROOT environment variable not set" -ForegroundColor Yellow
    Write-Host "Please specify vcpkg installation directory" -ForegroundColor Yellow
    Write-Host "Example: C:\vcpkg" -ForegroundColor Gray
    $VCPKG_ROOT = Read-Host "VCPKG_ROOT"
    
    if (-not (Test-Path $VCPKG_ROOT)) {
        Write-Host "`n[!] Error: vcpkg directory not found!" -ForegroundColor Red
        Write-Host "    Path: $VCPKG_ROOT" -ForegroundColor Red
        exit 1
    }
} else {
    $VCPKG_ROOT = $env:VCPKG_ROOT
    Write-Host "`nUsing VCPKG_ROOT: $VCPKG_ROOT" -ForegroundColor Green
}

$VCPKG_TOOLCHAIN = "$VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake"

if (-not (Test-Path $VCPKG_TOOLCHAIN)) {
    Write-Host "`n[!] Error: vcpkg toolchain file not found!" -ForegroundColor Red
    Write-Host "    Expected: $VCPKG_TOOLCHAIN" -ForegroundColor Red
    exit 1
}

# 检查ITK是否已安装
Write-Host "`n[*] Checking if ITK is installed via vcpkg..." -ForegroundColor Yellow
$vcpkgExe = "$VCPKG_ROOT\vcpkg.exe"
$installedPackages = & $vcpkgExe list | Select-String "itk"

if (-not $installedPackages) {
    Write-Host "[!] ITK not found in vcpkg" -ForegroundColor Red
    Write-Host "[*] Installing ITK... (this may take a while)" -ForegroundColor Yellow
    & $vcpkgExe install itk:x64-windows
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "`n[!] Failed to install ITK" -ForegroundColor Red
        exit 1
    }
} else {
    Write-Host "[✓] ITK is already installed" -ForegroundColor Green
}

# 检查是否存在build目录
if (Test-Path "build") {
    Write-Host "`n[?] Build directory exists. Remove it? (y/n)" -ForegroundColor Yellow
    $response = Read-Host
    if ($response -eq 'y' -or $response -eq 'Y') {
        Write-Host "[*] Removing old build directory..." -ForegroundColor Yellow
        Remove-Item -Recurse -Force "build"
    }
}

# 创建build目录
Write-Host "`n[1/2] Creating build directory..." -ForegroundColor Green
New-Item -ItemType Directory -Force -Path "build" | Out-Null
Set-Location "build"

# 配置CMake
Write-Host "`n[2/2] Configuring CMake with vcpkg toolchain..." -ForegroundColor Green
cmake -G "Visual Studio 16 2019" -A x64 `
    "-DCMAKE_TOOLCHAIN_FILE=$VCPKG_TOOLCHAIN" `
    ..

if ($LASTEXITCODE -ne 0) {
    Write-Host "`n[!] CMake configuration failed!" -ForegroundColor Red
    Set-Location ..
    exit 1
}

# 编译
Write-Host "`n[*] Building project (Release mode)..." -ForegroundColor Green
cmake --build . --config Release

if ($LASTEXITCODE -ne 0) {
    Write-Host "`n[!] Build failed!" -ForegroundColor Red
    Set-Location ..
    exit 1
}

# 完成
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "  Build Completed Successfully!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "`nExecutable location:" -ForegroundColor Yellow
Write-Host "  .\bin\Release\MIRegistration.exe" -ForegroundColor White
Write-Host "`nUsage:" -ForegroundColor Yellow
Write-Host "  .\bin\Release\MIRegistration.exe <fixed> <moving> <output>" -ForegroundColor White

Set-Location ..

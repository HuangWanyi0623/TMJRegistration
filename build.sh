#!/bin/bash
# 构建脚本 - Linux/macOS
# 使用方法: chmod +x build.sh && ./build.sh

echo "========================================"
echo "  MI Registration Build Script"
echo "========================================"

# 检查是否存在build目录
if [ -d "build" ]; then
    echo ""
    read -p "[?] Build directory exists. Remove it? (y/n): " response
    if [ "$response" = "y" ] || [ "$response" = "Y" ]; then
        echo "[*] Removing old build directory..."
        rm -rf build
    fi
fi

# 创建build目录
echo ""
echo "[1/3] Creating build directory..."
mkdir -p build
cd build

# 配置CMake
echo ""
echo "[2/3] Configuring CMake..."
echo "Please specify ITK installation directory"
echo "Example: /usr/local/lib/cmake/ITK-5.3"
read -p "ITK_DIR: " ITK_DIR

if [ ! -d "$ITK_DIR" ]; then
    echo ""
    echo "[!] Error: ITK directory not found!"
    echo "    Path: $ITK_DIR"
    cd ..
    exit 1
fi

cmake -DCMAKE_BUILD_TYPE=Release -DITK_DIR="$ITK_DIR" ..

if [ $? -ne 0 ]; then
    echo ""
    echo "[!] CMake configuration failed!"
    cd ..
    exit 1
fi

# 编译
echo ""
echo "[3/3] Building project..."
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)

if [ $? -ne 0 ]; then
    echo ""
    echo "[!] Build failed!"
    cd ..
    exit 1
fi

# 完成
echo ""
echo "========================================"
echo "  Build Completed Successfully!"
echo "========================================"
echo ""
echo "Executable location:"
echo "  ./bin/MIRegistration"
echo ""
echo "Usage:"
echo "  ./bin/MIRegistration <fixed> <moving> <output>"

cd ..

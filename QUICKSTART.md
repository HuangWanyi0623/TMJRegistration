# 快速开始指南

## 第一步: 安装vcpkg和ITK

### 1.1 安装vcpkg

```powershell
# 打开PowerShell,切换到合适的目录
cd C:\

# 克隆vcpkg仓库
git clone https://github.com/Microsoft/vcpkg.git

# 进入目录并运行bootstrap
cd vcpkg
.\bootstrap-vcpkg.bat

# 设置环境变量(临时)
$env:VCPKG_ROOT = "C:\vcpkg"

# 或永久设置(重启终端后生效)
[System.Environment]::SetEnvironmentVariable('VCPKG_ROOT', 'C:\vcpkg', 'User')
```

### 1.2 安装ITK

```powershell
# 在vcpkg目录下运行
.\vcpkg install itk:x64-windows

# 这个过程需要1-2小时,请耐心等待
# 安装完成后会看到类似输出:
# itk:x64-windows package is installed successfully
```

## 第二步: 编译项目

```powershell
# 返回项目目录
cd e:\图像处理\TMJ_12.6

# 在 PowerShell 中运行构建脚本
.\build.ps1

# 如果从 cmd 运行, 使用以下命令:
# powershell -ExecutionPolicy Bypass -File .\build.ps1

# 脚本会自动:
# 1. 检查vcpkg环境
# 2. 检查ITK是否已安装
# 3. 配置CMake
# 4. 编译项目
```

## 第三步: 准备测试数据

确保您有以下格式的医学图像:
- NRRD格式 (.nrrd)
- DICOM系列
- 或其他ITK支持格式

假设您有:
- `fixed_mri.nrrd` (固定图像)
- `moving_ct.nrrd` (移动图像)

## 第四步: 运行配准

### 方法1: 使用测试脚本

```powershell
.\test_registration.ps1
# 按提示输入文件路径
```

### 方法2: 直接运行

```powershell
.\build\bin\Release\MIRegistration.exe fixed_mri.nrrd moving_ct.nrrd output.nrrd
```

## 第五步: 查看结果

程序会输出:
1. 图像加载信息
2. 优化过程(实时)
3. 最终变换参数
4. 配准后图像保存位置

使用医学图像查看器(如3D Slicer, ITK-SNAP)查看配准结果。

---

## 常见问题

**Q: vcpkg安装太慢?**
- 可以使用代理或镜像
- 确保网络连接稳定

**Q: 编译失败?**
- 检查Visual Studio是否已安装
- 确保CMake版本 >= 3.16

**Q: 找不到图像文件?**
- 使用绝对路径
- 确保文件格式正确

**Q: 配准效果不好?**
- 检查图像是否已大致对齐
- 调整参数(在main.cpp中)

---

## 下一步

- 阅读 `README.md` 了解详细信息
- 修改 `src/main.cpp` 中的参数优化配准
- 根据需要改进算法

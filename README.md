# 基于互信息的3D医学图像配准工具

## 项目简介

这是一个用C++实现的基于Mattes互信息的3D医学图像配准工具,专门用于MRI和CT图像的精配准。本项目借鉴ITK的设计理念实现核心算法,便于后续改进和优化。

> **快速开始**: 如果您是第一次使用,请先阅读 [QUICKSTART.md](QUICKSTART.md)

### 主要特性

- **Mattes互信息度量**: 使用直方图方法估计联合概率分布
- **梯度下降优化器**: 自适应步长调整,自动收敛
- **刚体变换**: 支持3D平移和旋转 (6参数)
- **DICOM/NRRD支持**: 直接读取医学图像格式
- **高效采样**: 随机采样策略加速计算

## 环境要求

- Windows 10/11
- Visual Studio 2019或更高版本
- CMake 3.16+
- vcpkg (包管理器)

## 安装步骤

### 1. 安装vcpkg

```powershell
# 克隆vcpkg仓库
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg

# 运行bootstrap脚本
.\bootstrap-vcpkg.bat

# 设置环境变量
$env:VCPKG_ROOT = "C:\vcpkg"  # 替换为实际路径
```

### 2. 安装ITK

```powershell
# 使用vcpkg安装ITK (需要较长时间)
.\vcpkg install itk:x64-windows
```

### 3. 编译项目

```powershell
# 返回项目目录
cd e:\图像处理\TMJ_12.6

# 运行构建脚本
.\build.ps1
```

## 使用方法

### 基本用法

```powershell
.\build\bin\Release\MIRegistration.exe <固定图像> <移动图像> <输出图像>
```

### 示例

```powershell
# 使用DICOM系列或NRRD文件
.\build\bin\Release\MIRegistration.exe fixed_mri.nrrd moving_ct.nrrd output.nrrd
```

### 支持的图像格式

- **NRRD** (.nrrd, .nhdr) - 推荐格式
- **DICOM** (系列文件夹)
- **NIfTI** (.nii, .nii.gz)
- **MetaImage** (.mha, .mhd)

## 参数调整

在`main.cpp`中可修改配准参数:

```cpp
registration.SetNumberOfHistogramBins(50);        // 直方图bins: 30-100
registration.SetNumberOfSpatialSamples(10000);    // 采样点数: 5000-50000
registration.SetMaximumStepLength(1.0);           // 最大步长: 0.5-2.0
registration.SetMinimumStepLength(0.001);         // 最小步长: 0.0001-0.01
registration.SetNumberOfIterations(200);          // 最大迭代: 100-500
```

### 推荐配置

**快速配准** (1-2分钟):
- Bins: 30, Samples: 5000, Iterations: 100

**标准配准** (2-5分钟):
- Bins: 50, Samples: 10000, Iterations: 200

**高精度配准** (10-20分钟):
- Bins: 100, Samples: 50000, Iterations: 500

## 输出说明

程序会输出:
- 图像加载信息
- 优化过程 (每10次迭代)
- 最终变换参数 (平移mm, 旋转度)
- 配准耗时

## 常见问题

### 编译问题

**Q: vcpkg安装ITK很慢?**  
A: ITK包较大,首次安装需要1-2小时,请耐心等待。

**Q: 找不到vcpkg?**  
A: 确保设置了`VCPKG_ROOT`环境变量,或在运行`build.ps1`时手动输入路径。

### 运行问题

**Q: 配准不收敛?**  
A: 增加采样数到20000,增加迭代次数到300。

**Q: 速度太慢?**  
A: 减少采样数到5000,减少bins到30。

**Q: 精度不够?**  
A: 增加bins到80,增加采样数到20000。

## 后续改进方向

- [ ] 多分辨率策略
- [ ] 完整仿射变换 (12参数)
- [ ] OpenMP并行加速
- [ ] 更多优化器 (LBFGS, Powell)

## 参考资料

- ITK官网: https://itk.org/
- vcpkg: https://github.com/Microsoft/vcpkg
- Mattes等人论文: "Nonrigid multimodality image registration" (2001)

---

**版本**: 1.0.0  
**日期**: 2025-12-07  
**适用场景**: MRI-CT多模态配准, 刚体运动配准


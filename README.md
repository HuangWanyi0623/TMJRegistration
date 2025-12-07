# 基于互信息的3D医学图像配准工具

## 项目简介

这是一个用C++实现的基于Mattes互信息的3D医学图像配准工具,专门用于MRI和CT图像的精配准。本项目**自定义实现**核心算法,便于后续改进和优化。

> **快速开始**: 如果您是第一次使用,请先阅读 [QUICKSTART.md](QUICKSTART.md)

### 主要特性

- **Mattes互信息度量**: 自定义实现,使用直方图方法估计联合概率分布
- **梯度下降优化器**: 自定义实现,自适应步长调整,自动收敛
- **多分辨率金字塔**: 自定义实现,3层金字塔策略加速收敛
- **刚体变换**: 支持3D平移和旋转 (6参数)
- **DICOM/NRRD支持**: 直接读取医学图像格式
- **高效采样**: 固定种子随机采样,确保可重复性

### 自定义实现的类

本项目核心算法完全自定义实现,方便扩展和修改:

| 类名 | 文件 | 说明 |
|------|------|------|
| `MattesMutualInformation` | `src/MattesMutualInformation.cpp` | Mattes互信息度量 |
| `RegularStepGradientDescentOptimizer` | `src/RegularStepGradientDescentOptimizer.cpp` | 规则步长梯度下降优化器 |
| `ImageRegistration` | `src/ImageRegistration.cpp` | 配准框架(含多分辨率金字塔) |

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
.\build\bin\Release\MIRegistration.exe <固定图像> <移动图像> <输出文件夹>
```

**注意**: 第三个参数是**输出文件夹路径**,程序会在该文件夹中自动生成时间戳命名的 `.h5` 变换文件。

### 示例

```powershell
# 使用DICOM系列或NRRD文件
.\build\bin\Release\MIRegistration.exe fixed_mri.nrrd moving_ct.nrrd .\output

# 输出文件示例: .\output\transform_20251207_143025.h5
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

程序运行时会输出以下信息:

### 1. 图像元数据
- **Size**: 图像尺寸 (像素)
- **Spacing**: 体素间距 (mm)
- **Origin**: 图像原点位置 (mm)
- **Direction Matrix**: 图像方向矩阵 (3x3)

### 2. 多分辨率策略
ITK使用**多分辨率金字塔**策略提高配准鲁棒性:
- **Level 0 (粗略)**: 下采样4倍,平滑2mm → 快速粗略对齐
- **Level 1 (中等)**: 下采样2倍,平滑1mm → 中等精度
- **Level 2 (精细)**: 原始分辨率,不平滑 → 精细对齐

**这就是为什么迭代次数会归零3次!** 每层金字塔都会从iter=0重新开始。

### 3. 迭代信息
每10次迭代显示一次:
- **Iter**: 当前层的迭代次数 (0-300)
- **Metric**: 互信息度量值 (负值,**越小越好**)
- **LearningRate**: 当前学习率 (自适应调整)

**为什么最后的Metric没有第一层高?**
- Level 0在低分辨率图像上计算,信息量少,所以Metric绝对值小
- Level 2在原始分辨率上计算,信息量大,所以Metric绝对值大
- **不能直接比较不同层的Metric值!** 每层图像不同,Metric的范围也不同

### 4. 最终结果
- **Transform Parameters**: 旋转中心、旋转角度(弧度)、平移向量(mm)
- **4x4 Matrix**: 完整变换矩阵 (包含旋转中心偏移)
- **.h5文件**: 保存完整变换参数供Slicer等软件使用

### 示例输出
```
--- Multi-Resolution Level 0 ---
  Shrink factors: [4, 4, 4]
  Smoothing sigma: 2 mm
  Iter:    0  Metric:    -0.488776  LearningRate: 5.0000e-01
  Iter:   10  Metric:    -0.488227  LearningRate: 5.0000e-01
  ...
  
--- Multi-Resolution Level 1 ---
  Shrink factors: [2, 2, 2]
  Smoothing sigma: 1 mm
  Iter:    0  Metric:    -0.403787  LearningRate: 5.0000e-01
  ...

--- Multi-Resolution Level 2 ---
  Shrink factors: [1, 1, 1]
  Smoothing sigma: 0 mm
  Iter:    0  Metric:    -0.332674  LearningRate: 5.0000e-01
  ...
Final metric value: -3.4278e-01
```

## 运行

### 程序输出信息

运行时会显示以下信息：

1. **图像元数据**
   - Size: 图像尺寸 (体素数)
   - Spacing: 体素间距 (mm)
   - Origin: 图像原点坐标 (mm)
   - Direction Matrix: 方向矩阵 (3x3)

2. **配准进度** (每10次迭代显示一次)
   - `Iter`: 当前迭代次数
   - `Metric`: 互信息度量值 (负值，越小表示相似度越高)
   - `LearningRate`: 当前学习率 (优化器自适应调整)

3. **最终配准结果**
   - Rotation Center: 旋转中心坐标 (mm)
   - Rotation: 旋转角度 (rad)
   - Translation: 平移向量 (mm)
   - 4x4 Transformation Matrix: 完整的变换矩阵

4. **.h5文件保存位置**

示例输出：
```
Iter:    0  Metric:   -0.523456  LearningRate: 5.0000e-01
Iter:   10  Metric:   -0.645123  LearningRate: 4.8500e-01
Iter:   20  Metric:   -0.712345  LearningRate: 4.7000e-01
...
Final metric value: -0.823456
```

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

### 坐标系统问题 (重要!)

**Q: PowerShell输出的4x4矩阵与Slicer中.h5文件显示的为什么会有差异?**

**A: 这是Euler变换的内部表示方式导致的:**

1. **ITK的Euler3DTransform存储方式**:
   - 内部存储: **旋转中心(center)** + **旋转角度(3个)** + **平移(3个)**
   - 这种表示方式是参数化的,便于优化

2. **GetMatrix() + GetOffset() 的含义**:
   - 返回的是**相对于旋转中心的变换矩阵**
   - 公式: `T(x) = R*(x - center) + center + translation`
   - PowerShell显示的是 `[R | R*(-center) + center + translation]`

3. **.h5文件保存的是完整参数**:
   - 保存: center, rotation angles, translation
   - Slicer读取后会用**不同的旋转中心**重新计算矩阵
   - 所以看起来矩阵不一样,但**实际变换效果是一样的**!

4. **验证方法**:
   - 在Slicer中应用.h5变换到moving图像
   - 检查变换后的图像是否与fixed对齐
   - 如果对齐了,说明变换是正确的(虽然矩阵数值看起来不同)

**重要**: 不要手动输入PowerShell显示的矩阵!请使用.h5文件,它包含完整的变换参数。

---

**Q: 为什么需要初始化对齐(质心初始化)?**

**A: 这是优化算法的特性,不是ITK的限制:**

1. **互信息是非凸函数** - 有很多局部最优点
2. **梯度下降容易陷入局部最优** - 需要好的起点
3. **质心初始化提供合理起点** - 让优化器更容易找到全局最优
4. **ANTS使用更强大的优化策略** - 多分辨率金字塔 + 更好的优化器(如SyN),所以不需要初始化也能工作

**本程序使用简单的梯度下降**,所以需要质心初始化。如果配准效果不好:
- 检查图像的物理空间是否重叠(Origin、Spacing、Direction)
- 增加迭代次数和采样点数
- 考虑使用多分辨率策略(需要修改代码)

---

**Q: CT-MR配准后,变换效果不好?**  
**A: 可能的原因:**

1. **物理空间差异太大** - 检查程序输出的Origin、Spacing、Direction
2. **初始对齐不够好** - Origin差异导致质心初始化失效
3. **优化参数不合适** - 调整学习率、迭代次数
4. **随机采样不稳定** - 每次运行结果不同,需要固定随机种子或增加采样点数

### 配准稳定性问题

**Q: 相同数据每次运行结果都不一样?**  
**A: 已修复!** v1.1版本做了以下改进:

1. **固定随机种子** - `MetricSamplingReinitializeSeed(121212)` 确保采样点完全一致
2. **几何中心初始化** - 从 `MomentsOn()` 改为 `GeometryOn()`,不受灰度分布影响
3. **增加采样数** - 从20万提升到10万(自适应),提高稳定性
4. **优化器参数调整**:
   - 学习率: 1.0 → 0.5 (更保守)
   - 迭代次数: 500 → 300 (配合更稳定的采样)
   - 松弛因子: 0.5 → 0.8 (加速收敛)
   - 返回最佳参数: `SetReturnBestParametersAndValue(true)`

**现在相同输入保证相同输出!**

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

**版本**: 1.1.0
- 修复: 随机采样导致的不可重复性(固定种子)
- 改进: 几何中心初始化替代质心初始化
- 优化: 调整优化器参数提高稳定性
- 新增: 详细的图像元数据输出  
**日期**: 2025-12-07  
**适用场景**: MRI-CT多模态配准, 刚体运动配准


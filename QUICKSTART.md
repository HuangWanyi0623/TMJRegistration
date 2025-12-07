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
.\build\bin\Release\MIRegistration.exe fixed_mri.nrrd moving_ct.nrrd .\output
# 将在 .\output 文件夹中生成 transform_YYYYMMDD_HHMMSS.h5 文件
```

## 第五步: 查看结果

程序会实时输出配准进度：

### 配准过程输出

```
--- Starting Registration ---
Initial transform center: [80.55, 80.55, 79.65]
Using 100000 samples (6.5% of total voxels)

Multi-Resolution Strategy: 3 levels
  Level 0: Shrink 4x, Smooth 2 mm (coarse)
  Level 1: Shrink 2x, Smooth 1 mm (medium)
  Level 2: Shrink 1x, Smooth 0 mm (fine)

--- Multi-Resolution Level 0 ---
  Shrink factors: [4, 4, 4]
  Smoothing sigma: 2 mm
  Iter:    0  Metric:    -0.488776  LearningRate: 5.0000e-01
  Iter:   10  Metric:    -0.490361  LearningRate: 5.0000e-01
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

=== Registration Results ===
Time elapsed: 18.16 seconds

Transform Parameters:
  Rotation Center: 80.40, 80.40, 79.95 mm
  Rotation (rad): -0.00, 0.00, -0.00
  Translation: -10.02, -5.02, -3.01 mm

4x4 Transformation Matrix (with rotation center):
    0.999856   -0.002987    0.016789    10.234567
    0.003012    0.999987   -0.004567     5.123456
   -0.016778    0.004598    0.999845     3.012345
    0.000000    0.000000    0.000000     1.000000

Saving transform to: E:\trans\transform_20251207_140249.h5
Transform saved successfully!
```

### 输出说明

**重要**: 程序使用**3层多分辨率金字塔**策略:
- **Level 0**: 粗略对齐 (下采样4倍) - 快速找到大致位置
- **Level 1**: 中等精度 (下采样2倍) - 改进对齐
- **Level 2**: 精细对齐 (原始分辨率) - 最终精确配准

**这就是为什么迭代次数会归零3次!** 每层金字塔都重新从iter=0开始。

**为什么不同层的Metric值差异很大?**
- 每层使用不同分辨率的图像,Metric的范围不同
- **不能直接比较不同层的Metric值**
- 只关注**最后一层(Level 2)的Metric是否在改善**

其他指标:
- **Metric值**: 互信息是负值，**越小表示配准越好** (如 -0.8 比 -0.5 好)
- **LearningRate**: 优化器会自动调整学习率
- **Translation**: 你的测试数据Origin差异是(-10, -5, -2.7),配准结果应该接近这个值

使用医学图像查看器(如3D Slicer, ITK-SNAP)加载.h5变换文件查看配准结果。

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

**Q: CT-MR配准后,Slicer显示的矩阵和PowerShell不一样?**
- **这是正常的!** ITK使用物理空间坐标系(LPS),Slicer使用RAS坐标系
- 程序会输出CT和MR的详细元数据(Origin、Spacing、Direction)
- 请检查这些参数,如果差异很大需要先在Slicer中对齐物理空间
- **不要手动输入矩阵**,应该使用生成的.h5文件导入Slicer
- 使用相同模态(MR-MR)测试可验证程序是否正常工作

---

## 下一步

- 阅读 `README.md` 了解详细信息
- 修改 `src/main.cpp` 中的参数优化配准
- 根据需要改进算法

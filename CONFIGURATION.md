# 配准参数配置指南

本文档提供了不同场景下的推荐配置参数。

## 配置方式

### 方式1: JSON配置文件 (推荐) ⭐

使用JSON文件管理配准参数,支持命令行加载:

```powershell
.\build\bin\Release\MIRegistration.exe --config ".\config\Rigid.json" --initial "初始变换.h5" "fixed.nrrd" "moving.nrrd" "输出目录\"
```

**JSON配置示例** (`config/Rigid.json`):
```json
{
    "transformType": "Rigid",
    "numberOfHistogramBins": 32,
    "samplingPercentage": 0.25,
    "learningRate": [2.0, 1.0, 0.5, 0.1, 0.05],
    "minimumStepLength": 1e-6,
    "numberOfIterations": [1000, 500, 250, 100, 0],
    "relaxationFactor": 0.5,
    "gradientMagnitudeTolerance": 1e-6,
    "numberOfLevels": 5,
    "shrinkFactors": [12, 8, 4, 2, 1],
    "smoothingSigmas": [4.0, 3.0, 2.0, 1.0, 1.0],
    "useStratifiedSampling": true,
    "randomSeed": 121212
}
```

**关键特性**:
- ✅ **分层学习率**: `[2.0, 1.0, 0.5, 0.1, 0.05]` - 粗糙层大步长，精细层小步长
- ✅ 每层独立迭代次数: `[1000, 500, 250, 100, 0]`
- ✅ ANTs风格金字塔: 5层从粗到精
- ✅ 智能跳过策略: 最后一层设为0跳过全分辨率优化
- ✅ 采样百分比: 0.25 (25%采样)代替固定采样点数
- ✅ 更严格收敛: minimumStepLength=1e-6, gradientMagnitudeTolerance=1e-6

### 方式2: 代码中直接设置

在`main.cpp`中修改参数(传统方式):

```cpp
registration.SetNumberOfHistogramBins(50);
registration.SetNumberOfSpatialSamples(10000);
registration.SetMaximumStepLength(1.0);
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(200);
```

---

## 推荐配置 (基于实测)

### ⚡ ANTs风格5层金字塔 (最快,推荐) 

**适用场景**: 有粗配准初始化的精配准

```json
{
    "transformType": "Rigid",
    "numberOfHistogramBins": 32,
    "samplingPercentage": 0.25,
    "learningRate": [2.0, 1.0, 0.5, 0.1, 0.05],
    "minimumStepLength": 1e-6,
    "numberOfIterations": [1000, 500, 250, 100, 0],
    "numberOfLevels": 5,
    "shrinkFactors": [12, 8, 4, 2, 1],
    "smoothingSigmas": [4.0, 3.0, 2.0, 1.0, 1.0]
}
```

**性能**: 65-80秒 (153M体素, 8线程)  
**精度**: 高  
**优势**: 
- **分层学习率**: 粗糙层级(12x, 8x)使用大步长(2.0, 1.0)极快收敛
- 精细层级(4x, 2x)使用小步长(0.5, 0.1)避免过冲
- 跳过全分辨率优化节省大量时间
- 早期收敛自动停止,不浪费迭代
- 25%采样率提供更好的统计性
- 更严格的收敛阈值(1e-6)确保精度

---

## 传统配置示例 (代码方式)

### 1. 默认配置 (通用)

适用于大多数MRI-CT配准场景:

```cpp
registration.SetNumberOfHistogramBins(50);
registration.SetNumberOfSpatialSamples(10000);
registration.SetMaximumStepLength(1.0);
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(200);
```

**预期耗时**: 2-5分钟  
**适用场景**: 256³左右的标准医学图像  
**注意**: 不如ANTs 5层金字塔策略快

---

### 2. 快速配准 (低精度)

适用于预览或初步配准:

```cpp
registration.SetNumberOfHistogramBins(30);
registration.SetNumberOfSpatialSamples(5000);
registration.SetMaximumStepLength(2.0);
registration.SetMinimumStepLength(0.01);
registration.SetNumberOfIterations(100);
```

**预期耗时**: 30秒-1分钟  
**精度**: 中等,可能需要后续精配准

---

### 3. 高精度配准

适用于需要高精度的研究场景:

```cpp
registration.SetNumberOfHistogramBins(100);
registration.SetNumberOfSpatialSamples(50000);
registration.SetMaximumStepLength(0.5);
registration.SetMinimumStepLength(0.0001);
registration.SetNumberOfIterations(500);
```

**预期耗时**: 10-20分钟  
**精度**: 高,适合最终配准

---

### 4. 大图像配准 (512³+)

```cpp
registration.SetNumberOfHistogramBins(50);
registration.SetNumberOfSpatialSamples(100000);  // 增加采样
registration.SetMaximumStepLength(1.0);
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(300);
```

**预期耗时**: 15-30分钟  
**说明**: 大图像需要更多采样点保证覆盖

---

### 5. 小图像配准 (128³-)

```cpp
registration.SetNumberOfHistogramBins(30);
registration.SetNumberOfSpatialSamples(3000);    // 减少采样
registration.SetMaximumStepLength(1.0);
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(150);
```

**预期耗时**: 30秒-1分钟

---

### 6. 多模态配准 (MRI-CT)

```cpp
registration.SetNumberOfHistogramBins(80);       // 增加bins
registration.SetNumberOfSpatialSamples(20000);
registration.SetMaximumStepLength(0.5);          // 保守步长
registration.SetMinimumStepLength(0.0005);
registration.SetNumberOfIterations(300);
```

**说明**: 多模态图像对比度差异大,需要更细致的直方图

---

### 7. 同模态配准 (MRI-MRI, CT-CT)

```cpp
registration.SetNumberOfHistogramBins(30);       // 可以较少
registration.SetNumberOfSpatialSamples(10000);
registration.SetMaximumStepLength(1.5);          // 可以更激进
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(150);
```

**说明**: 同模态图像相似度高,收敛较快

---

## 参数说明

### NumberOfHistogramBins (直方图Bins数)
- **作用**: 控制概率密度估计的精细度
- **范围**: 20-200
- **建议**: 
  - 低对比度: 50-100
  - 高对比度: 20-50
- **影响**: 
  - 过小: 丢失信息
  - 过大: 计算慢,可能噪声敏感

### NumberOfSpatialSamples (空间采样点数)
- **作用**: 控制参与计算的体素数量
- **范围**: 1000-100000
- **建议**:
  - 快速预览: 3000-5000
  - 标准配准: 10000-20000
  - 高精度: 50000+
- **影响**:
  - 过小: 不稳定,可能不收敛
  - 过大: 计算慢

### MaximumStepLength (最大步长)
- **作用**: 初始搜索步长
- **范围**: 0.1-5.0
- **建议**:
  - 初始位置好: 0.5-1.0
  - 初始位置差: 1.0-2.0
  - 精细调整: 0.1-0.5
- **单位**: 与变换参数相关
  - 平移参数: mm
  - 旋转参数: 弧度

### MinimumStepLength (最小步长)
- **作用**: 收敛阈值
- **范围**: 0.00001-0.1
- **建议**:
  - 快速配准: 0.01
  - 标准配准: 0.001
  - 高精度: 0.0001

### NumberOfIterations (最大迭代次数)
- **作用**: 防止无限循环
- **范围**: 50-1000
- **建议**:
  - 快速配准: 50-100
  - 标准配准: 150-300
  - 高精度: 300-500
- **说明**: 通常会因其他收敛条件提前停止

---

## 调参建议流程

### 第一步: 使用默认参数测试

```cpp
// 使用默认配置运行一次
registration.SetNumberOfHistogramBins(50);
registration.SetNumberOfSpatialSamples(10000);
registration.SetMaximumStepLength(1.0);
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(200);
```

观察:
1. 是否收敛?
2. 收敛速度如何?
3. 最终结果精度?

### 第二步: 根据结果调整

#### 如果不收敛:
- 增加`NumberOfSpatialSamples` → 20000
- 增加`MaximumStepLength` → 2.0
- 增加`NumberOfIterations` → 300

#### 如果收敛太慢:
- 增加`MaximumStepLength` → 1.5
- 减少`NumberOfSpatialSamples` → 5000
- 调整`RelaxationFactor` → 0.6

#### 如果精度不够:
- 增加`NumberOfHistogramBins` → 80
- 增加`NumberOfSpatialSamples` → 20000
- 减小`MinimumStepLength` → 0.0005

#### 如果速度太慢:
- 减少`NumberOfSpatialSamples` → 5000
- 减少`NumberOfHistogramBins` → 30
- 减少`NumberOfIterations` → 100

### 第三步: 精细调优

记录每组参数的:
- 最终度量值
- 收敛迭代次数
- 总耗时
- 视觉评估结果

选择最佳平衡点。

---

## 代码模板

在`main.cpp`中修改参数的位置:

```cpp
// 创建配准对象
ImageRegistration registration;

// 设置图像
registration.SetFixedImage(fixedReader->GetOutput());
registration.SetMovingImage(movingReader->GetOutput());

// ============= 在这里修改参数 =============
registration.SetNumberOfHistogramBins(50);
registration.SetNumberOfSpatialSamples(10000);
registration.SetMaximumStepLength(1.0);
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(200);
// ========================================

// 执行配准
registration.Update();
```

---

## 故障排除速查表

| 症状 | 可能原因 | 解决方案 |
|-----|---------|---------|
| 不收敛 | 采样太少 | 增加samples到20000 |
| 收敛慢 | 步长太小 | 增加max_step到2.0 |
| 精度低 | Bins太少 | 增加bins到80 |
| 速度慢 | 采样太多 | 减少samples到5000 |
| 震荡 | 步长太大 | 减小max_step到0.5 |
| 局部极值 | 初始化差 | 手动初始化或粗配准 |

---

## JSON配置参数详解

### 基本参数

#### transformType
- **值**: "Rigid" 或 "Affine"
- **说明**: 变换类型
  - Rigid: 6参数刚体变换(旋转+平移)
  - Affine: 12参数仿射变换(旋转+平移+缩放+剪切)

#### numberOfHistogramBins
- **值**: 20-100 (整数)
- **推荐**: 32 (ANTs风格)
- **说明**: 互信息直方图bins数,影响概率密度估计精度

#### samplingPercentage
- **值**: 0.01-1.0 (小数)
- **推荐**: 0.25 (25%采样,精度与速度平衡)
- **说明**: 随机采样体素百分比,代替固定采样点数
- **优势**: 自动适应不同大小图像
- **平衡**: 
  - 0.1 (10%): 更快但可能不够稳定
  - 0.25 (25%): 推荐,统计性好
  - 0.5 (50%): 精度高但速度慢

### 优化器参数

#### learningRate
- **类型**: 
  - **单一值**: `0.1` (所有层使用相同学习率)
  - **分层数组**: `[2.0, 1.0, 0.5, 0.1, 0.05]` ⭐ **推荐**
- **推荐配置**:
  - ANTs 5层金字塔: `[2.0, 1.0, 0.5, 0.1, 0.05]`
  - 传统单一值: `0.1` (保守但可能较慢)
- **分层学习率原理**:
  - **粗糙层级** (12x, 8x): 高斯平滑后能量曲面接近凸函数，可以使用**大步长**（2.0, 1.0）快速逼近全局最优
  - **精细层级** (4x, 2x, 1x): 全分辨率下存在大量局部极小值，需要**小步长**（0.5, 0.1, 0.05）避免过冲
- **实验证据**:
  - 学习率 2.0 在粗糙层（Level 1-3）获得最高互信息值
  - 学习率 2.0 在精细层（Level 4）出现**过冲现象**：跳过最优点导致MI值反而降低
  - 分层策略兼顾速度（粗糙层快速收敛）和精度（精细层稳定优化）

**JSON配置示例**:
```json
{
    "learningRate": [2.0, 1.0, 0.5, 0.1, 0.05],  // 分层学习率
    "numberOfLevels": 5,
    "shrinkFactors": [12, 8, 4, 2, 1]
}
```

**向后兼容**: 仍支持单一值 `"learningRate": 0.1`，将自动应用到所有层级。

#### minimumStepLength
- **值**: 1e-8 到 1e-3
- **推荐**: 1e-6 (更严格,精度更高)
- **说明**: 收敛阈值,步长小于此值时停止

#### numberOfIterations
- **类型**: 整数数组 `[level0, level1, ...]`
- **推荐**: `[1000, 500, 250, 100, 0]` (ANTs 5层)
- **说明**: 每层金字塔的最大迭代次数
- **技巧**: 
  - 最后一层设为0可跳过全分辨率优化
  - 实际迭代次数通常因早期收敛而更少

#### relaxationFactor
- **值**: 0.3-0.9
- **推荐**: 0.5
- **说明**: 步长衰减因子,每次迭代步长 *= relaxationFactor

#### gradientMagnitudeTolerance
- **值**: 1e-8 到 1e-3
- **推荐**: 1e-6 (更严格)
- **说明**: 梯度大小收敛阈值

### 多分辨率参数

#### numberOfLevels
- **值**: 1-7 (整数)
- **推荐**: 5 (ANTs风格)
- **说明**: 金字塔层数

#### shrinkFactors
- **类型**: 整数数组 `[level0, level1, ...]`
- **推荐**: `[12, 8, 4, 2, 1]` (ANTs 5层)
- **说明**: 每层图像下采样因子
- **技巧**: 
  - 12x/8x粗糙层非常快(图像只有几千体素)
  - 1x = 全分辨率

#### smoothingSigmas
- **类型**: 浮点数组 `[level0, level1, ...]` (单位: mm)
- **推荐**: `[4.0, 3.0, 2.0, 1.0, 1.0]`
- **说明**: 每层图像高斯平滑标准差
- **作用**: 抑制噪声,扩大收敛域

### 采样参数

#### useStratifiedSampling
- **值**: true 或 false
- **推荐**: true
- **说明**: 使用分层采样确保空间分布均匀

#### randomSeed
- **值**: 任意整数
- **推荐**: 121212
- **说明**: 随机数种子,确保结果可重复

---

## ANTs策略深度解析

### 为什么5层金字塔这么快?

1. **极粗糙层级(12x, 8x)**
   - 12x下采样: 153M体素 → **88K体素** (1736倍缩减)
   - 8x下采样: 153M体素 → **298K体素** (513倍缩减)
   - **结果**: 这两层只需几秒钟

2. **早期收敛**
   - 优化器智能检测收敛
   - 粗糙层级通常20次迭代内收敛
   - 尽管配置了1000/500次,实际不会用完

3. **跳过全分辨率(1x)**
   - 设置 `numberOfIterations[4] = 0`
   - 只做初始化,不迭代优化
   - **节省**: 全分辨率优化通常需要100+秒

4. **渐进式优化**
   - 粗糙层快速找到大致位置
   - 中间层(4x, 2x)精细调整
   - 全分辨率层直接使用2x结果

### 何时使用ANTs策略?

✅ **适合**:
- 有粗配准初始化
- 需要极致速度
- 大图像(512³+)
- 多次配准任务

❌ **不适合**:
- 无初始化(从头配准)
- 需要最高精度(用传统3层+更多迭代)

### 如何调整ANTs参数?

**加快速度**:
- 增大粗糙层级: `[16, 12, 8, 4, 2]`
- 减少bins: 24 或 16
- 减少采样: 0.15 (15%)

**提高精度**:
- 保留1x层优化: `[1000, 500, 250, 100, 50]`
- 增加bins: 50
- 增加采样: 0.35 (35%)
- 更严格收敛: minimumStepLength=1e-7

---

**提示**: 建议保存每次实验的参数和结果,建立自己的参数数据库。

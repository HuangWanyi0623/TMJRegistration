# 算法实现技术文档

## 目录
1. [Mattes互信息算法详解](#1-mattes互信息算法详解)
2. [规则步长梯度下降优化器](#2-规则步长梯度下降优化器)
3. [刚体变换模型](#3-刚体变换模型)
4. [性能优化策略](#4-性能优化策略)
5. [参数调优指南](#5-参数调优指南)

---

## 1. Mattes互信息算法详解

### 1.1 理论基础

互信息(Mutual Information, MI)是一种基于信息论的相似性度量,定义为:

```
MI(A,B) = H(A) + H(B) - H(A,B)
       = Σ Σ p(a,b) log(p(a,b) / (p(a)·p(b)))
```

其中:
- `H(A)`, `H(B)` 是边缘熵
- `H(A,B)` 是联合熵
- `p(a,b)` 是联合概率
- `p(a)`, `p(b)` 是边缘概率

### 1.2 实现细节

#### 1.2.1 直方图构建

```cpp
// 初始化直方图
m_JointPDF[bins][bins] = 0
m_FixedMarginalPDF[bins] = 0
m_MovingMarginalPDF[bins] = 0

// 对每个采样点
for each sample_point:
    fixed_value = sample_point.value
    transformed_point = Transform(sample_point.position, parameters)
    moving_value = Interpolate(moving_image, transformed_point)
    
    fixed_bin = ComputeBin(fixed_value, fixed_range)
    moving_bin = ComputeBin(moving_value, moving_range)
    
    m_JointPDF[fixed_bin][moving_bin] += 1

// 归一化
for all bins:
    m_JointPDF[i][j] /= total_samples
    m_FixedMarginalPDF[i] = Σ m_JointPDF[i][j]
    m_MovingMarginalPDF[j] = Σ m_JointPDF[i][j]
```

#### 1.2.2 互信息计算

```cpp
MI = 0
for i in range(bins):
    for j in range(bins):
        if m_JointPDF[i][j] > epsilon:
            MI += m_JointPDF[i][j] * log(
                m_JointPDF[i][j] / 
                (m_FixedMarginalPDF[i] * m_MovingMarginalPDF[j])
            )
```

#### 1.2.3 梯度计算

使用有限差分法:

```cpp
for each parameter p[i]:
    p_plus = p
    p_plus[i] += delta
    
    gradient[i] = (GetValue(p_plus) - GetValue(p)) / delta
```

**优点**:
- 实现简单
- 数值稳定
- 易于调试

**缺点**:
- 需要多次计算度量值
- 计算成本较高

**改进方向**: 实现解析梯度(基于Parzen窗导数)

### 1.3 采样策略

#### 1.3.1 基础随机采样

**当前实现**: 固定图像随机均匀采样

```cpp
// 收集所有体素
all_voxels = GetAllVoxels(fixed_image)

// 随机打乱
Shuffle(all_voxels)

// 选择前N个
samples = all_voxels[0:N]
```

**优势**:
- 覆盖整个图像
- 减少计算量
- 避免偏差

#### 1.3.2 掩膜感知采样 (Mask-Aware Sampling)

**用途**: 局部配准 - 仅关注感兴趣区域(ROI)

**实现逻辑**:
```cpp
// 如果提供了掩膜
if (mask != null):
    // 1. 统计掩膜内体素数
    mask_voxels = CountVoxelsInMask(fixed_image, mask)
    
    // 2. 计算目标采样数 (百分比基于掩膜区域)
    num_samples = mask_voxels * sampling_percentage
    
    // 3. 采样时过滤掉掩膜外的点
    for each candidate_point:
        physical_point = IndexToPhysical(candidate_point)
        if mask.IsInside(physical_point):
            samples.add(candidate_point)
```

**关键修复 (2025-12-10/11)**:

**Bug #1: 采样数计算错误**
- **问题**: 采样数基于全图像计算,导致掩膜内过度采样
  - 示例: 全图1.5亿体素,掩膜覆盖15.8% (2400万),配置采样10%
  - 错误逻辑: 1.5亿 × 10% = 1534万样本
  - 实际效果: 1534万 / 2400万 = **63.3%过度采样**!
  
- **修复**: 采样数现在基于掩膜区域计算 (`MattesMutualInformation.cpp` 第89-122行)
  - 正确逻辑: 2400万 × 10% = 240万样本
  - 采样率: 240万 / 2400万 = **10%合理采样**
  
**Bug #2: 显示信息错误**
- **问题**: 输出仍显示 `Using 15341231 samples (10.0% of total voxels)` (基于全图)
- **修复**: 现在正确显示 `Using 2423974 samples (10.0% of mask region, 24239738 voxels in mask)`
  - 修复位置: `ImageRegistration.cpp` 第1089-1130行
  - 优化: 复用 `LoadFixedMask()` 统计的体素数,避免重复遍历图像

**实现要点**:
1. `LoadFixedMask()`: 加载时统计并保存 `m_MaskVoxelCount`
2. `MattesMutualInformation::Initialize()`: 使用Mask体素数计算采样数
3. `ImageRegistration::Update()`: 打印采样信息时使用保存的Mask体素数

**医学应用场景**:
- TMJ (颞下颌关节) 局部配准: CBCT全头 + MRI局部切片
- 脊椎局部配准: 忽略周围软组织
- 肿瘤局部配准: 仅关注病灶区域

**性能影响**:
- 计算量降低: 仅处理ROI内采样点
- 精度提升: 避免无关区域干扰
- 收敛速度: 取决于ROI复杂度

**可选策略**:
- 重要性采样(高梯度区域)
- 分层采样(保证各区域都有采样)
- 多掩膜加权采样

---

## 2. 规则步长梯度下降优化器

### 2.1 算法流程

```
初始化:
    position = initial_position
    step_length = maximum_step_length

for iteration in range(max_iterations):
    # 计算梯度
    gradient = ComputeGradient(position)
    gradient_magnitude = ||gradient||
    
    # 检查收敛
    if gradient_magnitude < tolerance:
        break
    
    # 归一化梯度
    gradient = gradient / gradient_magnitude
    
    # 尝试更新
    new_position = position - step_length * gradient
    new_value = ComputeValue(new_position)
    
    # 判断是否接受
    if new_value < current_value:
        position = new_position
        current_value = new_value
    else:
        # 减小步长
        step_length *= relaxation_factor
    
    # 检查最小步长
    if step_length < minimum_step_length:
        break
```

### 2.2 关键参数

| 参数 | 默认值 | 推荐范围 | 作用 |
|-----|--------|---------|-----|
| maximum_step_length | 1.0 | 0.1-2.0 | 初始搜索步长 |
| minimum_step_length | 0.001 | 0.0001-0.01 | 收敛阈值 |
| relaxation_factor | 0.5 | 0.3-0.7 | 步长缩减率 |
| gradient_tolerance | 1e-4 | 1e-5-1e-3 | 梯度幅值阈值 |

### 2.3 收敛判据

程序在满足以下任一条件时停止:

1. **步长过小**: `step_length < minimum_step_length`
2. **梯度过小**: `||gradient|| < gradient_tolerance`
3. **达到最大迭代次数**: `iteration >= max_iterations`

### 2.4 改进方向

1. **自适应步长**: 根据度量值变化自动调整
2. **动量法**: 添加历史梯度信息
3. **Line Search**: 沿梯度方向最优化步长
4. **更高级优化器**:
   - LBFGS (拟牛顿法)
   - Powell (无梯度方法)
   - Adam (自适应矩估计)

---

## 3. 刚体变换模型

### 3.1 参数定义

当前实现使用6参数刚体变换:

```cpp
parameters[0] = tx  // X轴平移 (mm)
parameters[1] = ty  // Y轴平移 (mm)
parameters[2] = tz  // Z轴平移 (mm)
parameters[3] = rx  // 绕X轴旋转 (弧度)
parameters[4] = ry  // 绕Y轴旋转 (弧度)
parameters[5] = rz  // 绕Z轴旋转 (弧度)
```

### 3.2 变换公式

使用ZYX欧拉角表示旋转:

```cpp
// 旋转矩阵
Rx = [1    0      0   ]
     [0  cos(rx) -sin(rx)]
     [0  sin(rx)  cos(rx)]

Ry = [cos(ry)  0  sin(ry)]
     [0        1    0     ]
     [-sin(ry) 0  cos(ry)]

Rz = [cos(rz) -sin(rz)  0]
     [sin(rz)  cos(rz)  0]
     [0        0        1]

// 组合旋转矩阵
R = Rz * Ry * Rx

// 完整变换
[x']   [R00 R01 R02]   [x]   [tx]
[y'] = [R10 R11 R12] * [y] + [ty]
[z']   [R20 R21 R22]   [z]   [tz]
```

### 3.3 扩展为仿射变换

要实现12参数仿射变换,需要修改参数结构:

```cpp
parameters[0-8]  = 旋转+缩放+剪切矩阵的9个元素
parameters[9-11] = 平移向量的3个元素
```

### 3.4 约束与正则化

对于医学图像配准,通常需要添加约束:

1. **旋转角度限制**: |rx|, |ry|, |rz| < 30°
2. **平移距离限制**: |tx|, |ty|, |tz| < 50mm
3. **正则化项**: 防止过度形变

---

## 4. 性能优化策略

### 4.1 当前性能特征

- **时间复杂度**: O(iterations × samples × bins)
- **空间复杂度**: O(bins² + samples)
- **典型耗时**: 
  - 图像加载: 1-5秒
  - 初始化: 1-2秒
  - 每次迭代: 0.5-2秒
  - 总耗时: 2-5分钟 (200次迭代)

### 4.2 已实现的优化

1. **随机采样**: 减少计算点数
2. **缓存采样点**: 避免重复采样
3. **高效插值**: 使用线性插值
4. **内存预分配**: 避免动态分配

### 4.3 可实现的优化

#### 4.3.1 多线程并行化

```cpp
#pragma omp parallel for
for (int i = 0; i < m_SamplePoints.size(); ++i)
{
    // 并行计算每个采样点的贡献
}
```

**预期加速**: 4-8倍 (取决于CPU核心数)

#### 4.3.2 SIMD向量化

```cpp
// 使用SSE/AVX指令加速
__m256 values = _mm256_load_ps(data);
__m256 result = _mm256_add_ps(values, offset);
```

**预期加速**: 2-4倍

#### 4.3.3 GPU加速

使用CUDA实现核心计算:

- 并行采样
- 并行直方图构建
- 并行插值

**预期加速**: 10-50倍

#### 4.3.4 多分辨率策略

```
Level 3 (低分辨率): 快速粗配准
    ↓ (上采样参数)
Level 2 (中分辨率): 精细调整
    ↓ (上采样参数)
Level 1 (原始分辨率): 最终优化
```

**优势**:
- 避免局部极值
- 加快收敛速度
- 提高配准精度

---

## 5. 参数调优指南

### 5.1 根据图像特征调整

#### 大图像 (>256³)
```cpp
SetNumberOfSpatialSamples(50000);     // 增加采样
SetNumberOfHistogramBins(100);        // 增加bins
SetNumberOfIterations(300);           // 增加迭代
```

#### 小图像 (<128³)
```cpp
SetNumberOfSpatialSamples(5000);      // 减少采样
SetNumberOfHistogramBins(30);         // 减少bins
SetNumberOfIterations(100);           // 减少迭代
```

### 5.2 根据图像对比度调整

#### 高对比度 (CT-CT, MRI-MRI)
```cpp
SetNumberOfHistogramBins(30);         // 较少bins即可
SetMaximumStepLength(2.0);            // 可以更激进
```

#### 低对比度/多模态 (MRI-CT)
```cpp
SetNumberOfHistogramBins(50-100);     // 需要更多bins
SetMaximumStepLength(0.5);            // 更保守的步长
```

### 5.3 调试技巧

#### 观察收敛曲线

在优化器中添加观察者:

```cpp
optimizer->SetObserver([](int iter, double value, const auto& params) {
    // 保存每次迭代的度量值
    // 绘制收敛曲线
});
```

**期望曲线**: 单调递减,最终趋于平稳

**异常情况**:
- 震荡: 步长过大
- 平台: 陷入局部极值
- 发散: 初始位置过差

#### 可视化中间结果

每N次迭代保存配准图像:

```cpp
if (iteration % 20 == 0) {
    SaveIntermediateResult(iteration);
}
```

### 5.4 常见问题解决

| 问题 | 可能原因 | 解决方案 |
|-----|---------|---------|
| 配准不收敛 | 初始位置差异大 | 先做粗配准或手动初始化 |
| 结果不准确 | 采样点太少 | 增加采样数到20000+ |
| 速度太慢 | 采样点太多 | 减少到5000-10000 |
| 陷入局部极值 | 初始步长太小 | 增大maximum_step_length |
| 结果不稳定 | 步长减小太快 | 调整relaxation_factor |

---

## 附录: ITK源码对照

### A.1 核心类对应关系

| 本项目 | ITK源码 |
|-------|---------|
| MattesMutualInformation | itkMattesMutualInformationImageToImageMetricv4 |
| RegularStepGradientDescentOptimizer | itkRegularStepGradientDescentOptimizerv4 |
| ImageRegistration | itkImageRegistrationMethodv4 |

### A.2 主要差异

1. **度量计算**: ITK使用B样条Parzen窗,本实现使用简化的直方图方法
2. **梯度计算**: ITK支持解析梯度,本实现使用有限差分
3. **变换模型**: ITK支持多种变换,本实现目前仅支持刚体变换
4. **多分辨率**: ITK内置多分辨率框架,本实现暂未包含

### A.3 算法等效性

在以下条件下,本实现与ITK等效:
- 使用相同的采样点
- 使用相同的直方图bins数
- 使用相同的优化参数
- 使用相同的插值方法

---

**文档版本**: 1.0  
**最后更新**: 2025-12-07

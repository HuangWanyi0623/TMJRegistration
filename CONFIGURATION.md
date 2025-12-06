# 配准参数配置示例

本文档提供了不同场景下的推荐配置参数。

## 1. 默认配置 (通用)

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

---

## 2. 快速配准 (低精度)

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

## 3. 高精度配准

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

## 4. 大图像配准 (512³+)

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

## 5. 小图像配准 (128³-)

```cpp
registration.SetNumberOfHistogramBins(30);
registration.SetNumberOfSpatialSamples(3000);    // 减少采样
registration.SetMaximumStepLength(1.0);
registration.SetMinimumStepLength(0.001);
registration.SetNumberOfIterations(150);
```

**预期耗时**: 30秒-1分钟

---

## 6. 多模态配准 (MRI-CT)

```cpp
registration.SetNumberOfHistogramBins(80);       // 增加bins
registration.SetNumberOfSpatialSamples(20000);
registration.SetMaximumStepLength(0.5);          // 保守步长
registration.SetMinimumStepLength(0.0005);
registration.SetNumberOfIterations(300);
```

**说明**: 多模态图像对比度差异大,需要更细致的直方图

---

## 7. 同模态配准 (MRI-MRI, CT-CT)

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

**提示**: 建议保存每次实验的参数和结果,建立自己的参数数据库。

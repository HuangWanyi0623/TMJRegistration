#include "MattesMutualInformation.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkGradientImageFilter.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numeric>

// ============================================================================
// 构造函数和析构函数
// ============================================================================

MattesMutualInformation::MattesMutualInformation()
    : m_NumberOfHistogramBins(50)
    , m_NumberOfSpatialSamples(0) // will be computed from sampling percentage by default
    , m_SamplingPercentage(0.10)   // default 10%
    , m_RandomSeed(121212)
    , m_UseFixedSeed(true)
    , m_UseStratifiedSampling(true)  // 默认使用分层采样
    , m_NumberOfValidSamples(0)
    , m_NumberOfParameters(6)  // 默认刚体6参数
    , m_FixedImageMin(0.0)
    , m_FixedImageMax(1.0)
    , m_MovingImageMin(0.0)
    , m_MovingImageMax(1.0)
    , m_FixedImageBinSize(1.0)
    , m_MovingImageBinSize(1.0)
    , m_CurrentValue(0.0)
    , m_NumberOfThreads(std::thread::hardware_concurrency())  // 自动检测CPU核心数
{
    m_Interpolator = InterpolatorType::New();
    m_RandomGenerator.seed(m_RandomSeed);
    
    // 初始化梯度插值器
    for (int i = 0; i < 3; ++i)
    {
        m_GradientInterpolators[i] = InterpolatorType::New();
    }
    
    if (m_NumberOfThreads == 0) m_NumberOfThreads = 1;
    
    std::cout << "[Metric] Multi-threading enabled: " << m_NumberOfThreads << " threads" << std::endl;
}

MattesMutualInformation::~MattesMutualInformation()
{
}

// ============================================================================
// 设置图像
// ============================================================================

void MattesMutualInformation::SetFixedImage(ImageType::Pointer fixedImage)
{
    m_FixedImage = fixedImage;
}

void MattesMutualInformation::SetMovingImage(ImageType::Pointer movingImage)
{
    m_MovingImage = movingImage;
    m_Interpolator->SetInputImage(m_MovingImage);
}

// ============================================================================
// 初始化
// ============================================================================

void MattesMutualInformation::Initialize()
{
    if (!m_FixedImage || !m_MovingImage)
    {
        throw std::runtime_error("Fixed or moving image not set");
    }

    // 使用固定种子确保可重复性
    if (m_UseFixedSeed)
    {
        m_RandomGenerator.seed(m_RandomSeed);
    }

    // 计算图像强度范围
    ComputeImageExtrema();
    
    // 计算bin大小
    m_FixedImageBinSize = (m_FixedImageMax - m_FixedImageMin) / m_NumberOfHistogramBins;
    m_MovingImageBinSize = (m_MovingImageMax - m_MovingImageMin) / m_NumberOfHistogramBins;

    // 基于百分比计算采样数量 (考虑掩膜区域)
    if (m_SamplingPercentage > 0.0 && m_NumberOfSpatialSamples == 0)
    {
        auto region = m_FixedImage->GetLargestPossibleRegion();
        auto size = region.GetSize();
        unsigned long totalVoxels = static_cast<unsigned long>(size[0]) * size[1] * size[2];
        
        // 如果有掩膜,需要先计算掩膜内的体素数
        if (m_FixedImageMask.IsNotNull())
        {
            unsigned long maskVoxels = 0;
            using IteratorType = itk::ImageRegionConstIteratorWithIndex<ImageType>;
            IteratorType it(m_FixedImage, region);
            for (it.GoToBegin(); !it.IsAtEnd(); ++it)
            {
                ImageType::PointType physicalPoint;
                m_FixedImage->TransformIndexToPhysicalPoint(it.GetIndex(), physicalPoint);
                if (m_FixedImageMask->IsInsideInWorldSpace(physicalPoint))
                {
                    ++maskVoxels;
                }
            }
            // 使用掩膜内体素数计算采样数
            m_NumberOfSpatialSamples = static_cast<unsigned int>(maskVoxels * m_SamplingPercentage + 0.5);
            if (m_Verbose)
            {
                std::cout << "[Metric] Mask-aware sampling: " << maskVoxels << " voxels in mask, "
                          << m_NumberOfSpatialSamples << " samples (" 
                          << (m_SamplingPercentage * 100.0) << "% of mask)" << std::endl;
            }
        }
        else
        {
            // 没有掩膜时使用全图像体素数
            m_NumberOfSpatialSamples = static_cast<unsigned int>(totalVoxels * m_SamplingPercentage + 0.5);
        }
    }

    if (m_Verbose)
    {
        std::cout << "[Metric Debug] Initialize: samplingPercentage=" << m_SamplingPercentage
                  << " numberOfSpatialSamples=" << m_NumberOfSpatialSamples << std::endl;
    }

    // 在固定图像上采样
    SampleFixedImage();

    // 计算移动图像梯度供解析梯度使用
    ComputeMovingImageGradient();
    if (m_Verbose)
    {
        std::cout << "[Metric Debug] Computed moving image gradient" << std::endl;
    }

    // 初始化直方图
    m_JointPDF.resize(m_NumberOfHistogramBins, 
                     std::vector<double>(m_NumberOfHistogramBins, 0.0));
    m_FixedImageMarginalPDF.resize(m_NumberOfHistogramBins, 0.0);
    m_MovingImageMarginalPDF.resize(m_NumberOfHistogramBins, 0.0);
    
    // 初始化梯度PDF存储 (根据参数数量动态分配)
    m_JointPDFDerivatives.resize(m_NumberOfParameters);
    for (auto& paramDerivative : m_JointPDFDerivatives)
    {
        paramDerivative.resize(m_NumberOfHistogramBins,
                              std::vector<double>(m_NumberOfHistogramBins, 0.0));
    }
}

void MattesMutualInformation::ReinitializeSampling()
{
    if (m_UseFixedSeed)
    {
        m_RandomGenerator.seed(m_RandomSeed);
    }
    SampleFixedImage();
}

// ============================================================================
// 图像极值计算
// ============================================================================

void MattesMutualInformation::ComputeImageExtrema()
{
    using IteratorType = itk::ImageRegionConstIteratorWithIndex<ImageType>;
    
    // 计算固定图像的最小最大值
    IteratorType it(m_FixedImage, m_FixedImage->GetLargestPossibleRegion());
    it.GoToBegin();
    m_FixedImageMin = it.Get();
    m_FixedImageMax = it.Get();
    
    while (!it.IsAtEnd())
    {
        double value = it.Get();
        if (value < m_FixedImageMin) m_FixedImageMin = value;
        if (value > m_FixedImageMax) m_FixedImageMax = value;
        ++it;
    }

    // 计算移动图像的最小最大值
    IteratorType it2(m_MovingImage, m_MovingImage->GetLargestPossibleRegion());
    it2.GoToBegin();
    m_MovingImageMin = it2.Get();
    m_MovingImageMax = it2.Get();
    
    while (!it2.IsAtEnd())
    {
        double value = it2.Get();
        if (value < m_MovingImageMin) m_MovingImageMin = value;
        if (value > m_MovingImageMax) m_MovingImageMax = value;
        ++it2;
    }
    
    // 稍微扩展范围避免边界问题
    double fixedRange = m_FixedImageMax - m_FixedImageMin;
    double movingRange = m_MovingImageMax - m_MovingImageMin;
    m_FixedImageMin -= fixedRange * 0.001;
    m_FixedImageMax += fixedRange * 0.001;
    m_MovingImageMin -= movingRange * 0.001;
    m_MovingImageMax += movingRange * 0.001;
}

// ============================================================================
// 计算移动图像梯度 (用于解析梯度)
// ============================================================================

void MattesMutualInformation::ComputeMovingImageGradient()
{
    // 使用ITK的GradientImageFilter计算图像梯度
    using GradientFilterType = itk::GradientImageFilter<ImageType, float, float>;
    using GradientOutputType = GradientFilterType::OutputImageType;
    
    auto gradientFilter = GradientFilterType::New();
    gradientFilter->SetInput(m_MovingImage);
    gradientFilter->SetUseImageSpacing(true);  // 考虑物理spacing
    gradientFilter->Update();
    
    GradientOutputType::Pointer gradientImage = gradientFilter->GetOutput();
    
    // 将梯度向量图像分解为3个标量图像(便于插值)
    ImageType::RegionType region = m_MovingImage->GetLargestPossibleRegion();
    
    for (int dim = 0; dim < 3; ++dim)
    {
        m_MovingImageGradient[dim] = ImageType::New();
        m_MovingImageGradient[dim]->SetRegions(region);
        m_MovingImageGradient[dim]->SetSpacing(m_MovingImage->GetSpacing());
        m_MovingImageGradient[dim]->SetOrigin(m_MovingImage->GetOrigin());
        m_MovingImageGradient[dim]->SetDirection(m_MovingImage->GetDirection());
        m_MovingImageGradient[dim]->Allocate();
    }
    
    // 复制梯度分量
    using GradientIteratorType = itk::ImageRegionConstIteratorWithIndex<GradientOutputType>;
    GradientIteratorType gradIt(gradientImage, region);
    
    for (gradIt.GoToBegin(); !gradIt.IsAtEnd(); ++gradIt)
    {
        auto index = gradIt.GetIndex();
        auto gradientVector = gradIt.Get();
        
        for (int dim = 0; dim < 3; ++dim)
        {
            m_MovingImageGradient[dim]->SetPixel(index, gradientVector[dim]);
        }
    }
    
    // 设置梯度插值器
    for (int dim = 0; dim < 3; ++dim)
    {
        m_GradientInterpolators[dim]->SetInputImage(m_MovingImageGradient[dim]);
    }
}

// ============================================================================
// 采样策略
// ============================================================================

void MattesMutualInformation::SampleFixedImage()
{
    // 随机采样固定图像
    using IteratorType = itk::ImageRegionConstIteratorWithIndex<ImageType>;
    ImageType::RegionType region = m_FixedImage->GetLargestPossibleRegion();

    // 收集所有有效点的索引
    std::vector<ImageType::IndexType> allIndices;
    IteratorType it(m_FixedImage, region);

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        // 如果设置了掩膜,只收集掩膜内部的索引
        if (m_FixedImageMask.IsNotNull())
        {
            ImageType::PointType physicalPoint;
            m_FixedImage->TransformIndexToPhysicalPoint(it.GetIndex(), physicalPoint);
            if (!m_FixedImageMask->IsInsideInWorldSpace(physicalPoint))
            {
                continue;  // 跳过掩膜外的点
            }
        }
        allIndices.push_back(it.GetIndex());
    }

    // 使用固定种子随机打乱
    std::shuffle(allIndices.begin(), allIndices.end(), m_RandomGenerator);

    unsigned int numSamples = std::min(m_NumberOfSpatialSamples, 
                                      static_cast<unsigned int>(allIndices.size()));
    
    m_SamplePoints.clear();
    m_SamplePoints.reserve(numSamples);

    for (unsigned int i = 0; i < numSamples; ++i)
    {
        SamplePoint sample;
        m_FixedImage->TransformIndexToPhysicalPoint(allIndices[i], sample.fixedPoint);
        sample.fixedValue = m_FixedImage->GetPixel(allIndices[i]);
        // compute parzen window index and weights for fixed image intensity
        double fixedContinuousIndex = ComputeFixedImageContinuousIndex(sample.fixedValue);
        int fixedStartIndex;
        ComputeBSplineWeights(fixedContinuousIndex, fixedStartIndex, sample.fixedBSplineWeights);
        sample.fixedParzenWindowIndex = fixedStartIndex;
        m_SamplePoints.push_back(sample);
    }

    m_NumberOfValidSamples = static_cast<unsigned int>(m_SamplePoints.size());

    if (m_Verbose)
    {
        std::cout << "[Metric Debug] SampleFixedImage: numSamples=" << numSamples
                  << " validSamples=" << m_NumberOfValidSamples 
                  << " maskEnabled=" << (m_FixedImageMask.IsNotNull() ? "Yes" : "No") << std::endl;
    }
}

void MattesMutualInformation::SampleFixedImageStratified()
{
    // 分层均匀采样: 将图像划分为网格,在每个格子中随机采样
    // 这确保采样点在空间上均匀分布
    
    ImageType::RegionType region = m_FixedImage->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    
    // 计算网格划分
    // 目标: 采样数 ≈ gridX * gridY * gridZ * samplesPerCell
    double totalVoxels = static_cast<double>(size[0] * size[1] * size[2]);
    double samplingRatio = static_cast<double>(m_NumberOfSpatialSamples) / totalVoxels;
    
    // 每个维度的网格数 (立方根近似)
    unsigned int gridDivisions = static_cast<unsigned int>(
        std::cbrt(static_cast<double>(m_NumberOfSpatialSamples) / 8.0));
    gridDivisions = std::max(gridDivisions, 2u);
    gridDivisions = std::min(gridDivisions, 50u);
    
    unsigned int gridX = std::min(gridDivisions, static_cast<unsigned int>(size[0]));
    unsigned int gridY = std::min(gridDivisions, static_cast<unsigned int>(size[1]));
    unsigned int gridZ = std::min(gridDivisions, static_cast<unsigned int>(size[2]));
    
    // 每个格子的采样数
    unsigned int totalCells = gridX * gridY * gridZ;
    unsigned int samplesPerCell = std::max(1u, m_NumberOfSpatialSamples / totalCells);
    
    // 每个格子的尺寸
    unsigned int cellSizeX = size[0] / gridX;
    unsigned int cellSizeY = size[1] / gridY;
    unsigned int cellSizeZ = size[2] / gridZ;
    
    m_SamplePoints.clear();
    m_SamplePoints.reserve(m_NumberOfSpatialSamples);
    
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (unsigned int gz = 0; gz < gridZ; ++gz)
    {
        for (unsigned int gy = 0; gy < gridY; ++gy)
        {
            for (unsigned int gx = 0; gx < gridX; ++gx)
            {
                // 当前格子的起始索引
                unsigned int startX = gx * cellSizeX;
                unsigned int startY = gy * cellSizeY;
                unsigned int startZ = gz * cellSizeZ;
                
                // 当前格子的结束索引
                unsigned int endX = (gx == gridX - 1) ? size[0] : (gx + 1) * cellSizeX;
                unsigned int endY = (gy == gridY - 1) ? size[1] : (gy + 1) * cellSizeY;
                unsigned int endZ = (gz == gridZ - 1) ? size[2] : (gz + 1) * cellSizeZ;
                
                // 在格子内随机采样
                for (unsigned int s = 0; s < samplesPerCell; ++s)
                {
                    if (m_SamplePoints.size() >= m_NumberOfSpatialSamples)
                        break;
                    
                    // 随机索引
                    ImageType::IndexType index;
                    index[0] = startX + static_cast<unsigned int>(dist(m_RandomGenerator) * (endX - startX));
                    index[1] = startY + static_cast<unsigned int>(dist(m_RandomGenerator) * (endY - startY));
                    index[2] = startZ + static_cast<unsigned int>(dist(m_RandomGenerator) * (endZ - startZ));
                    
                    // 边界检查
                    index[0] = std::min(index[0], static_cast<ImageType::IndexType::IndexValueType>(size[0] - 1));
                    index[1] = std::min(index[1], static_cast<ImageType::IndexType::IndexValueType>(size[1] - 1));
                    index[2] = std::min(index[2], static_cast<ImageType::IndexType::IndexValueType>(size[2] - 1));
                    
                    // 转换为物理坐标
                    ImageType::PointType physicalPoint;
                    m_FixedImage->TransformIndexToPhysicalPoint(index, physicalPoint);
                    
                    // 掩膜检查: 如果设置了掩膜,只采样掩膜内部的点
                    if (m_FixedImageMask.IsNotNull())
                    {
                        if (!m_FixedImageMask->IsInsideInWorldSpace(physicalPoint))
                        {
                            continue;  // 跳过掩膜外的点
                        }
                    }
                    
                    SamplePoint sample;
                    sample.fixedPoint = physicalPoint;
                    sample.fixedValue = m_FixedImage->GetPixel(index);
                    
                    // 预计算固定图像B样条权重
                    double fixedContinuousIndex = ComputeFixedImageContinuousIndex(sample.fixedValue);
                    ComputeBSplineWeights(fixedContinuousIndex, 
                                         sample.fixedParzenWindowIndex,
                                         sample.fixedBSplineWeights);
                    
                    m_SamplePoints.push_back(sample);
                }
            }
        }
    }
}

void MattesMutualInformation::SampleFixedImageRandom()
{
    // 原始随机采样方法(作为备选)
    using IteratorType = itk::ImageRegionConstIteratorWithIndex<ImageType>;
    ImageType::RegionType region = m_FixedImage->GetLargestPossibleRegion();
    
    std::vector<ImageType::IndexType> allIndices;
    IteratorType it(m_FixedImage, region);
    
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        // 如果设置了掩膜,只收集掩膜内部的索引
        if (m_FixedImageMask.IsNotNull())
        {
            ImageType::PointType physicalPoint;
            m_FixedImage->TransformIndexToPhysicalPoint(it.GetIndex(), physicalPoint);
            if (!m_FixedImageMask->IsInsideInWorldSpace(physicalPoint))
            {
                continue;  // 跳过掩膜外的点
            }
        }
        allIndices.push_back(it.GetIndex());
    }

    std::shuffle(allIndices.begin(), allIndices.end(), m_RandomGenerator);

    unsigned int numSamples = std::min(m_NumberOfSpatialSamples, 
                                      static_cast<unsigned int>(allIndices.size()));
    
    m_SamplePoints.clear();
    m_SamplePoints.reserve(numSamples);

    for (unsigned int i = 0; i < numSamples; ++i)
    {
        SamplePoint sample;
        m_FixedImage->TransformIndexToPhysicalPoint(allIndices[i], sample.fixedPoint);
        sample.fixedValue = m_FixedImage->GetPixel(allIndices[i]);
        
        // 预计算固定图像B样条权重
        double fixedContinuousIndex = ComputeFixedImageContinuousIndex(sample.fixedValue);
        ComputeBSplineWeights(fixedContinuousIndex, 
                             sample.fixedParzenWindowIndex,
                             sample.fixedBSplineWeights);
        
        m_SamplePoints.push_back(sample);
    }
}

// ============================================================================
// B样条计算
// ============================================================================

double MattesMutualInformation::EvaluateCubicBSpline(double u) const
{
    // 三次B样条基函数 (ITK使用的标准形式)
    // 定义在 [-2, 2] 区间
    double absU = std::abs(u);
    
    if (absU < 1.0)
    {
        // |u| < 1: (4 - 6u^2 + 3|u|^3) / 6
        return (4.0 - 6.0 * absU * absU + 3.0 * absU * absU * absU) / 6.0;
    }
    else if (absU < 2.0)
    {
        // 1 <= |u| < 2: (2 - |u|)^3 / 6
        double tmp = 2.0 - absU;
        return (tmp * tmp * tmp) / 6.0;
    }
    else
    {
        return 0.0;
    }
}

double MattesMutualInformation::EvaluateCubicBSplineDerivative(double u) const
{
    // 三次B样条的导数
    double absU = std::abs(u);
    double sign = (u >= 0) ? 1.0 : -1.0;
    
    if (absU < 1.0)
    {
        // d/du [(4 - 6u^2 + 3|u|^3) / 6] = (-12u + 9u|u|) / 6 = u(-2 + 1.5|u|)
        return sign * (-2.0 * absU + 1.5 * absU * absU);
    }
    else if (absU < 2.0)
    {
        // d/du [(2 - |u|)^3 / 6] = -sign * (2 - |u|)^2 / 2
        double tmp = 2.0 - absU;
        return -sign * (tmp * tmp) / 2.0;
    }
    else
    {
        return 0.0;
    }
}

void MattesMutualInformation::ComputeBSplineWeights(
    double continuousIndex, 
    int& startIndex, 
    std::array<double, 4>& weights) const
{
    // 计算B样条的起始索引和4个权重
    // B样条在 [startIndex-1, startIndex+2] 范围内有值
    
    startIndex = static_cast<int>(std::floor(continuousIndex)) - 1;
    
    for (int i = 0; i < 4; ++i)
    {
        double u = continuousIndex - (startIndex + i);
        weights[i] = EvaluateCubicBSpline(u);
    }
}

void MattesMutualInformation::ComputeBSplineDerivativeWeights(
    double continuousIndex, 
    int& startIndex, 
    std::array<double, 4>& derivativeWeights) const
{
    startIndex = static_cast<int>(std::floor(continuousIndex)) - 1;
    
    for (int i = 0; i < 4; ++i)
    {
        double u = continuousIndex - (startIndex + i);
        derivativeWeights[i] = EvaluateCubicBSplineDerivative(u);
    }
}

// ============================================================================
// 强度到连续索引转换
// ============================================================================

double MattesMutualInformation::ComputeFixedImageContinuousIndex(double value) const
{
    // 将图像强度值转换为直方图的连续索引
    // 留出2个bin的padding
    const double padding = 2.0;
    return padding + (value - m_FixedImageMin) / m_FixedImageBinSize;
}

double MattesMutualInformation::ComputeMovingImageContinuousIndex(double value) const
{
    const double padding = 2.0;
    return padding + (value - m_MovingImageMin) / m_MovingImageBinSize;
}

// ============================================================================
// 核心计算: 联合PDF和导数
// ============================================================================

void MattesMutualInformation::ComputeJointPDFAndDerivatives()
{
    if (!m_Transform)
    {
        throw std::runtime_error("Transform not set in metric");
    }
    
    if (!m_JacobianFunction)
    {
        throw std::runtime_error("Jacobian function not set in metric");
    }

    // 清空直方图
    for (auto& row : m_JointPDF)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    std::fill(m_FixedImageMarginalPDF.begin(), m_FixedImageMarginalPDF.end(), 0.0);
    std::fill(m_MovingImageMarginalPDF.begin(), m_MovingImageMarginalPDF.end(), 0.0);
    
    // 清空梯度PDF
    for (auto& paramDerivative : m_JointPDFDerivatives)
    {
        for (auto& row : paramDerivative)
        {
            std::fill(row.begin(), row.end(), 0.0);
        }
    }

    m_NumberOfValidSamples = 0;
    std::vector<std::array<double, 3>> jacobian;

    // 遍历所有采样点
    for (const auto& sample : m_SamplePoints)
    {
        // 使用变换将固定图像点变换到移动图像空间
        ImageType::PointType transformedPoint = m_Transform->TransformPoint(sample.fixedPoint);

        // 检查变换后的点是否在移动图像范围内
        if (!m_Interpolator->IsInsideBuffer(transformedPoint))
        {
            continue;
        }

        // 插值获取移动图像值
        double movingValue = m_Interpolator->Evaluate(transformedPoint);
        
        // 计算移动图像的连续索引和B样条权重
        double movingContinuousIndex = ComputeMovingImageContinuousIndex(movingValue);
        int movingStartIndex;
        std::array<double, 4> movingBSplineWeights;
        ComputeBSplineWeights(movingContinuousIndex, movingStartIndex, movingBSplineWeights);
        
        // 计算移动图像B样条导数权重
        std::array<double, 4> movingBSplineDerivativeWeights;
        int tempStartIndex;
        ComputeBSplineDerivativeWeights(movingContinuousIndex, tempStartIndex, movingBSplineDerivativeWeights);
        
        // 获取移动图像梯度
        std::array<double, 3> movingGradient = {0.0, 0.0, 0.0};
        for (int dim = 0; dim < 3; ++dim)
        {
            if (m_GradientInterpolators[dim]->IsInsideBuffer(transformedPoint))
            {
                movingGradient[dim] = m_GradientInterpolators[dim]->Evaluate(transformedPoint);
            }
        }
        
        // 使用外部提供的雅可比函数计算变换雅可比矩阵
        m_JacobianFunction(sample.fixedPoint, jacobian);
        
        // 计算 dm/dp = gradient_M^T * dT/dp
        // dm/dp[k] = sum_d (gradient_M[d] * jacobian[k][d])
        std::vector<double> dmDp(m_NumberOfParameters, 0.0);
        for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
        {
            dmDp[k] = 0.0;
            for (int d = 0; d < 3; ++d)
            {
                dmDp[k] += movingGradient[d] * jacobian[k][d];
            }
            // 转换为bin索引的导数
            dmDp[k] /= m_MovingImageBinSize;
        }

        // Print debug info for first few samples
        if (m_Verbose && m_NumberOfValidSamples < 5)
        {
            std::cout << "[Metric Debug] Sample " << m_NumberOfValidSamples << " fixedVal=" << sample.fixedValue
                      << " movedVal=" << movingValue << std::endl;
            std::cout << "  movingGradient: " << movingGradient[0] << " " << movingGradient[1] << " " << movingGradient[2] << std::endl;
            std::cout << "  jacobian[0]: " << jacobian[0][0] << " " << jacobian[0][1] << " " << jacobian[0][2] << std::endl;
            std::cout << "  dmDp (first 6): ";
            for (unsigned int k = 0; k < std::min<unsigned int>(6, dmDp.size()); ++k)
            {
                std::cout << dmDp[k] << " ";
            }
            std::cout << std::endl;
        }
        
        // 累加到联合PDF和导数PDF
        for (int fi = 0; fi < 4; ++fi)
        {
            int fixedBin = sample.fixedParzenWindowIndex + fi;
            if (fixedBin < 0 || fixedBin >= static_cast<int>(m_NumberOfHistogramBins))
                continue;
                
            double fixedWeight = sample.fixedBSplineWeights[fi];
            
            for (int mi = 0; mi < 4; ++mi)
            {
                int movingBin = movingStartIndex + mi;
                if (movingBin < 0 || movingBin >= static_cast<int>(m_NumberOfHistogramBins))
                    continue;
                
                double movingWeight = movingBSplineWeights[mi];
                double movingDerivWeight = movingBSplineDerivativeWeights[mi];
                
                // 联合PDF贡献
                double jointContribution = fixedWeight * movingWeight;
                m_JointPDF[fixedBin][movingBin] += jointContribution;
                
                // 导数PDF贡献
                // dP/dp = B_fixed * dB_moving/dm * dm/dp
                for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
                {
                    double derivContribution = fixedWeight * movingDerivWeight * dmDp[k];
                    m_JointPDFDerivatives[k][fixedBin][movingBin] += derivContribution;
                }
            }
        }
        
        m_NumberOfValidSamples++;
    }

    // 归一化为概率分布
    if (m_NumberOfValidSamples > 0)
    {
        double normFactor = 1.0 / static_cast<double>(m_NumberOfValidSamples);
        
        for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
        {
            for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
            {
                m_JointPDF[i][j] *= normFactor;
                m_FixedImageMarginalPDF[i] += m_JointPDF[i][j];
                m_MovingImageMarginalPDF[j] += m_JointPDF[i][j];
                
                // 归一化导数
                for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
                {
                    m_JointPDFDerivatives[k][i][j] *= normFactor;
                }
            }
        }
    }

    if (m_Verbose)
    {
        unsigned int nonZeroBins = 0;
        for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
        {
            for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
            {
                if (m_JointPDF[i][j] > 0.0)
                    ++nonZeroBins;
            }
        }
        double fillRatio = static_cast<double>(nonZeroBins) / (m_NumberOfHistogramBins * m_NumberOfHistogramBins);
        std::cout << "[Metric Debug] JointPDF filled bins: " << nonZeroBins << " (" << fillRatio << ")" << std::endl;
    }
}

// ============================================================================
// 多线程版本: 核心计算 - 联合PDF和导数 (高性能)
// ============================================================================

void MattesMutualInformation::ComputePDFRange(
    size_t startIdx, 
    size_t endIdx, 
    ThreadLocalHistograms& localHist)
{
    std::vector<std::array<double, 3>> jacobian;
    
    // 处理分配给此线程的采样点
    for (size_t sampleIdx = startIdx; sampleIdx < endIdx; ++sampleIdx)
    {
        const auto& sample = m_SamplePoints[sampleIdx];
        
        // 使用变换将固定图像点变换到移动图像空间
        ImageType::PointType transformedPoint = m_Transform->TransformPoint(sample.fixedPoint);

        // 检查变换后的点是否在移动图像范围内
        if (!m_Interpolator->IsInsideBuffer(transformedPoint))
        {
            continue;
        }

        // 插值获取移动图像值
        double movingValue = m_Interpolator->Evaluate(transformedPoint);
        
        // 计算移动图像的连续索引和B样条权重
        double movingContinuousIndex = ComputeMovingImageContinuousIndex(movingValue);
        int movingStartIndex;
        std::array<double, 4> movingBSplineWeights;
        ComputeBSplineWeights(movingContinuousIndex, movingStartIndex, movingBSplineWeights);
        
        // 计算移动图像B样条导数权重
        std::array<double, 4> movingBSplineDerivativeWeights;
        int tempStartIndex;
        ComputeBSplineDerivativeWeights(movingContinuousIndex, tempStartIndex, movingBSplineDerivativeWeights);
        
        // 获取移动图像梯度
        std::array<double, 3> movingGradient = {0.0, 0.0, 0.0};
        for (int dim = 0; dim < 3; ++dim)
        {
            if (m_GradientInterpolators[dim]->IsInsideBuffer(transformedPoint))
            {
                movingGradient[dim] = m_GradientInterpolators[dim]->Evaluate(transformedPoint);
            }
        }
        
        // 使用外部提供的雅可比函数计算变换雅可比矩阵
        m_JacobianFunction(sample.fixedPoint, jacobian);
        
        // 计算 dm/dp = gradient_M^T * dT/dp
        std::vector<double> dmDp(m_NumberOfParameters, 0.0);
        for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
        {
            dmDp[k] = 0.0;
            for (int d = 0; d < 3; ++d)
            {
                dmDp[k] += movingGradient[d] * jacobian[k][d];
            }
            dmDp[k] /= m_MovingImageBinSize;
        }
        
        // 累加到局部线程的联合PDF和导数PDF
        for (int fi = 0; fi < 4; ++fi)
        {
            int fixedBin = sample.fixedParzenWindowIndex + fi;
            if (fixedBin < 0 || fixedBin >= static_cast<int>(m_NumberOfHistogramBins))
                continue;
                
            double fixedWeight = sample.fixedBSplineWeights[fi];
            
            for (int mi = 0; mi < 4; ++mi)
            {
                int movingBin = movingStartIndex + mi;
                if (movingBin < 0 || movingBin >= static_cast<int>(m_NumberOfHistogramBins))
                    continue;
                
                double movingWeight = movingBSplineWeights[mi];
                double movingDerivWeight = movingBSplineDerivativeWeights[mi];
                
                // 联合PDF贡献
                double jointContribution = fixedWeight * movingWeight;
                localHist.jointPDF[fixedBin][movingBin] += jointContribution;
                
                // 导数PDF贡献
                for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
                {
                    double derivContribution = fixedWeight * movingDerivWeight * dmDp[k];
                    localHist.jointPDFDerivatives[k][fixedBin][movingBin] += derivContribution;
                }
            }
        }
        
        localHist.validSamples++;
    }
}

void MattesMutualInformation::ComputeJointPDFAndDerivativesThreaded()
{
    if (!m_Transform)
    {
        throw std::runtime_error("Transform not set in metric");
    }
    
    if (!m_JacobianFunction)
    {
        throw std::runtime_error("Jacobian function not set in metric");
    }

    // 清空全局直方图
    for (auto& row : m_JointPDF)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    std::fill(m_FixedImageMarginalPDF.begin(), m_FixedImageMarginalPDF.end(), 0.0);
    std::fill(m_MovingImageMarginalPDF.begin(), m_MovingImageMarginalPDF.end(), 0.0);
    
    for (auto& paramDerivative : m_JointPDFDerivatives)
    {
        for (auto& row : paramDerivative)
        {
            std::fill(row.begin(), row.end(), 0.0);
        }
    }

    // 创建线程局部直方图
    std::vector<ThreadLocalHistograms> threadHistograms;
    threadHistograms.reserve(m_NumberOfThreads);
    for (unsigned int t = 0; t < m_NumberOfThreads; ++t)
    {
        threadHistograms.emplace_back(m_NumberOfHistogramBins, m_NumberOfParameters);
    }

    // 将采样点分配给各个线程
    size_t totalSamples = m_SamplePoints.size();
    size_t samplesPerThread = totalSamples / m_NumberOfThreads;
    
    std::vector<std::thread> threads;
    threads.reserve(m_NumberOfThreads);
    
    for (unsigned int t = 0; t < m_NumberOfThreads; ++t)
    {
        size_t startIdx = t * samplesPerThread;
        size_t endIdx = (t == m_NumberOfThreads - 1) ? totalSamples : (t + 1) * samplesPerThread;
        
        threads.emplace_back([this, startIdx, endIdx, &threadHistograms, t]() {
            this->ComputePDFRange(startIdx, endIdx, threadHistograms[t]);
        });
    }
    
    // 等待所有线程完成
    for (auto& thread : threads)
    {
        thread.join();
    }
    
    // 合并所有线程的局部直方图到全局直方图 (Reduce阶段)
    m_NumberOfValidSamples = 0;
    for (const auto& localHist : threadHistograms)
    {
        m_NumberOfValidSamples += localHist.validSamples;
        
        for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
        {
            for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
            {
                m_JointPDF[i][j] += localHist.jointPDF[i][j];
                
                for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
                {
                    m_JointPDFDerivatives[k][i][j] += localHist.jointPDFDerivatives[k][i][j];
                }
            }
        }
    }

    // 归一化为概率分布
    if (m_NumberOfValidSamples > 0)
    {
        double normFactor = 1.0 / static_cast<double>(m_NumberOfValidSamples);
        
        for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
        {
            for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
            {
                m_JointPDF[i][j] *= normFactor;
                m_FixedImageMarginalPDF[i] += m_JointPDF[i][j];
                m_MovingImageMarginalPDF[j] += m_JointPDF[i][j];
                
                for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
                {
                    m_JointPDFDerivatives[k][i][j] *= normFactor;
                }
            }
        }
    }

    if (m_Verbose)
    {
        unsigned int nonZeroBins = 0;
        for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
        {
            for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
            {
                if (m_JointPDF[i][j] > 0.0)
                    ++nonZeroBins;
            }
        }
        double fillRatio = static_cast<double>(nonZeroBins) / (m_NumberOfHistogramBins * m_NumberOfHistogramBins);
        std::cout << "[Metric Debug - Multithreaded] JointPDF filled bins: " << nonZeroBins 
                  << " (" << fillRatio << "), Valid samples: " << m_NumberOfValidSamples << std::endl;
    }
}

// ============================================================================
// 计算互信息值
// ============================================================================

double MattesMutualInformation::ComputeMutualInformation()
{
    double mutualInformation = 0.0;
    const double epsilon = 1e-16;

    for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
    {
        double fixedProb = m_FixedImageMarginalPDF[i];
        
        if (fixedProb < epsilon)
            continue;

        for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
        {
            double movingProb = m_MovingImageMarginalPDF[j];
            double jointProb = m_JointPDF[i][j];

            if (jointProb < epsilon || movingProb < epsilon)
                continue;

            // MI = sum( P(f,m) * log( P(f,m) / (P(f) * P(m)) ) )
            mutualInformation += jointProb * std::log(jointProb / (fixedProb * movingProb));
        }
    }

    return mutualInformation;
}

// ============================================================================
// 计算解析梯度
// ============================================================================

void MattesMutualInformation::ComputeAnalyticalGradient(ParametersType& derivative)
{
    // 解析梯度公式:
    // dMI/dp = sum_f sum_m [ dP(f,m)/dp * (1 + log(P(f,m) / P(m))) ]
    //        = sum_f sum_m [ dP(f,m)/dp * log(P(f,m) / P(m)) ]
    //          + sum_f sum_m [ dP(f,m)/dp ]
    // 
    // 第二项 = d(sum P)/dp = 0 (因为概率和为1)
    // 所以:
    // dMI/dp = sum_f sum_m [ dP(f,m)/dp * log(P(f,m) / P(m)) ]
    
    const double epsilon = 1e-16;
    derivative.resize(m_NumberOfParameters, 0.0);
    
    for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
    {
        for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
        {
            double jointProb = m_JointPDF[i][j];
            double movingProb = m_MovingImageMarginalPDF[j];
            
            if (jointProb < epsilon || movingProb < epsilon)
                continue;
            
            // log(P(f,m) / P(m))
            double logTerm = std::log(jointProb / movingProb);
            
            for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
            {
                derivative[k] += m_JointPDFDerivatives[k][i][j] * logTerm;
            }
        }
    }
    
    // 因为我们最小化负互信息,所以梯度取负
    for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
    {
        derivative[k] = -derivative[k];
    }

    if (m_Verbose)
    {
        double maxAbs = 0.0, sumAbs = 0.0;
        for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
        {
            double a = std::abs(derivative[k]);
            if (a > maxAbs) maxAbs = a;
            sumAbs += a;
        }
        double meanAbs = sumAbs / m_NumberOfParameters;
        std::cout << "[Metric Debug] Gradient stats: maxAbs=" << maxAbs << " meanAbs=" << meanAbs << std::endl;
        std::cout << "  gradient: ";
        for (unsigned int k = 0; k < m_NumberOfParameters; ++k)
            std::cout << std::fixed << std::setprecision(6) << derivative[k] << " ";
        std::cout << std::endl;
    }
}

// ============================================================================
// 公共接口
// ============================================================================

double MattesMutualInformation::GetValue()
{
    ComputeJointPDFAndDerivativesThreaded();
    double mi = ComputeMutualInformation();
    m_CurrentValue = -mi;  // 返回负值,因为我们要最小化
    return m_CurrentValue;
}

void MattesMutualInformation::GetDerivative(ParametersType& derivative)
{
    // 先计算PDF(如果还没计算的话)
    ComputeJointPDFAndDerivativesThreaded();
    ComputeAnalyticalGradient(derivative);
}

void MattesMutualInformation::GetValueAndDerivative(double& value, ParametersType& derivative)
{
    ComputeJointPDFAndDerivativesThreaded();
    value = -ComputeMutualInformation();
    m_CurrentValue = value;
    ComputeAnalyticalGradient(derivative);
}

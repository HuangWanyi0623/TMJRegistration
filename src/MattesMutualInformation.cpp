#include "MattesMutualInformation.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

MattesMutualInformation::MattesMutualInformation()
    : m_NumberOfHistogramBins(50)
    , m_NumberOfSpatialSamples(100000)
    , m_RandomSeed(121212)
    , m_UseFixedSeed(true)
    , m_FixedImageMin(0.0)
    , m_FixedImageMax(1.0)
    , m_MovingImageMin(0.0)
    , m_MovingImageMax(1.0)
    , m_CurrentValue(0.0)
{
    m_Interpolator = InterpolatorType::New();
    m_RandomGenerator.seed(m_RandomSeed);
}

MattesMutualInformation::~MattesMutualInformation()
{
}

void MattesMutualInformation::SetFixedImage(ImageType::Pointer fixedImage)
{
    m_FixedImage = fixedImage;
}

void MattesMutualInformation::SetMovingImage(ImageType::Pointer movingImage)
{
    m_MovingImage = movingImage;
    m_Interpolator->SetInputImage(m_MovingImage);
}

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

    // 在固定图像上采样
    SampleFixedImage();

    // 初始化直方图
    m_JointPDF.resize(m_NumberOfHistogramBins, 
                     std::vector<double>(m_NumberOfHistogramBins, 0.0));
    m_FixedImageMarginalPDF.resize(m_NumberOfHistogramBins, 0.0);
    m_MovingImageMarginalPDF.resize(m_NumberOfHistogramBins, 0.0);
}

void MattesMutualInformation::ReinitializeSampling()
{
    // 重新使用固定种子
    if (m_UseFixedSeed)
    {
        m_RandomGenerator.seed(m_RandomSeed);
    }
    SampleFixedImage();
}

void MattesMutualInformation::ComputeImageExtrema()
{
    // 计算固定图像的最小最大值
    using IteratorType = itk::ImageRegionConstIteratorWithIndex<ImageType>;
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
}

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
        m_FixedImage->TransformIndexToPhysicalPoint(allIndices[i], sample.point);
        sample.fixedValue = m_FixedImage->GetPixel(allIndices[i]);
        m_SamplePoints.push_back(sample);
    }
}

void MattesMutualInformation::ComputeJointPDF()
{
    if (!m_Transform)
    {
        throw std::runtime_error("Transform not set in metric");
    }

    // 清空直方图
    for (auto& row : m_JointPDF)
    {
        std::fill(row.begin(), row.end(), 0.0);
    }
    std::fill(m_FixedImageMarginalPDF.begin(), m_FixedImageMarginalPDF.end(), 0.0);
    std::fill(m_MovingImageMarginalPDF.begin(), m_MovingImageMarginalPDF.end(), 0.0);

    unsigned int validSamples = 0;

    // 遍历所有采样点
    for (const auto& sample : m_SamplePoints)
    {
        // 使用变换将固定图像点变换到移动图像空间
        ImageType::PointType transformedPoint = m_Transform->TransformPoint(sample.point);

        // 检查变换后的点是否在移动图像范围内
        if (!m_Interpolator->IsInsideBuffer(transformedPoint))
        {
            continue;
        }

        // 插值获取移动图像值
        double movingValue = m_Interpolator->Evaluate(transformedPoint);

        // 计算直方图bin索引
        int fixedBin = ComputeFixedImageBin(sample.fixedValue);
        int movingBin = ComputeMovingImageBin(movingValue);

        if (fixedBin >= 0 && fixedBin < static_cast<int>(m_NumberOfHistogramBins) &&
            movingBin >= 0 && movingBin < static_cast<int>(m_NumberOfHistogramBins))
        {
            m_JointPDF[fixedBin][movingBin] += 1.0;
            validSamples++;
        }
    }

    // 归一化为概率分布
    if (validSamples > 0)
    {
        double normalizationFactor = 1.0 / validSamples;
        
        for (unsigned int i = 0; i < m_NumberOfHistogramBins; ++i)
        {
            for (unsigned int j = 0; j < m_NumberOfHistogramBins; ++j)
            {
                m_JointPDF[i][j] *= normalizationFactor;
                m_FixedImageMarginalPDF[i] += m_JointPDF[i][j];
                m_MovingImageMarginalPDF[j] += m_JointPDF[i][j];
            }
        }
    }
}

double MattesMutualInformation::ComputeMutualInformation()
{
    double mutualInformation = 0.0;
    const double epsilon = 1e-12; // 避免log(0)

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

            // MI = sum( P(i,j) * log( P(i,j) / (P(i) * P(j)) ) )
            mutualInformation += jointProb * std::log(jointProb / (fixedProb * movingProb));
        }
    }

    return mutualInformation;
}

double MattesMutualInformation::GetValue()
{
    ComputeJointPDF();
    double mi = ComputeMutualInformation();
    m_CurrentValue = -mi; // 返回负值,因为我们要最小化
    return m_CurrentValue;
}

void MattesMutualInformation::GetDerivative(ParametersType& derivative)
{
    // 使用有限差分法计算梯度
    // 对于Euler3DTransform: 参数0-2是旋转(rad), 参数3-5是平移(mm)
    const unsigned int numParams = 6;
    derivative.resize(numParams);
    
    // 不同参数类型使用不同的扰动量
    // 旋转参数用小扰动(弧度), 平移参数用较大扰动(毫米)
    const double rotationDelta = 1e-4;    // ~0.006 degrees
    const double translationDelta = 0.01;  // 0.01 mm

    // 保存当前变换参数
    auto currentParams = m_Transform->GetParameters();
    
    // 使用中心差分法,更准确
    for (unsigned int i = 0; i < numParams; ++i)
    {
        double delta = (i < 3) ? rotationDelta : translationDelta;
        
        // 正向扰动
        auto forwardParams = currentParams;
        forwardParams[i] += delta;
        m_Transform->SetParameters(forwardParams);
        double forwardValue = GetValue();
        
        // 负向扰动
        auto backwardParams = currentParams;
        backwardParams[i] -= delta;
        m_Transform->SetParameters(backwardParams);
        double backwardValue = GetValue();
        
        // 中心差分
        derivative[i] = (forwardValue - backwardValue) / (2.0 * delta);
    }

    // 恢复原始参数
    m_Transform->SetParameters(currentParams);
    // 重新计算当前值
    m_CurrentValue = GetValue();
}

void MattesMutualInformation::GetValueAndDerivative(double& value, ParametersType& derivative)
{
    value = GetValue();
    GetDerivative(derivative);
}

int MattesMutualInformation::ComputeFixedImageBin(double value) const
{
    if (m_FixedImageMax <= m_FixedImageMin)
        return 0;
        
    if (value <= m_FixedImageMin)
        return 0;
    if (value >= m_FixedImageMax)
        return m_NumberOfHistogramBins - 1;

    double range = m_FixedImageMax - m_FixedImageMin;
    int bin = static_cast<int>((value - m_FixedImageMin) / range * m_NumberOfHistogramBins);
    
    return std::min(bin, static_cast<int>(m_NumberOfHistogramBins - 1));
}

int MattesMutualInformation::ComputeMovingImageBin(double value) const
{
    if (m_MovingImageMax <= m_MovingImageMin)
        return 0;
        
    if (value <= m_MovingImageMin)
        return 0;
    if (value >= m_MovingImageMax)
        return m_NumberOfHistogramBins - 1;

    double range = m_MovingImageMax - m_MovingImageMin;
    int bin = static_cast<int>((value - m_MovingImageMin) / range * m_NumberOfHistogramBins);
    
    return std::min(bin, static_cast<int>(m_NumberOfHistogramBins - 1));
}

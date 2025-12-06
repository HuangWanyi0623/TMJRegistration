#include "MattesMutualInformation.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>

MattesMutualInformation::MattesMutualInformation()
    : m_NumberOfHistogramBins(50)
    , m_NumberOfSpatialSamples(10000)
    , m_FixedImageMin(0.0)
    , m_FixedImageMax(1.0)
    , m_MovingImageMin(0.0)
    , m_MovingImageMax(1.0)
{
    m_Interpolator = InterpolatorType::New();
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

    // 计算图像强度范围
    ComputeImageExtrema();

    // 在固定图像上采样
    SampleFixedImage();

    // 初始化直方图
    m_JointPDF.resize(m_NumberOfHistogramBins, 
                     std::vector<double>(m_NumberOfHistogramBins, 0.0));
    m_FixedImageMarginalPDF.resize(m_NumberOfHistogramBins, 0.0);
    m_MovingImageMarginalPDF.resize(m_NumberOfHistogramBins, 0.0);

    std::cout << "MattesMutualInformation initialized:" << std::endl;
    std::cout << "  Number of histogram bins: " << m_NumberOfHistogramBins << std::endl;
    std::cout << "  Number of samples: " << m_SamplePoints.size() << std::endl;
    std::cout << "  Fixed image range: [" << m_FixedImageMin << ", " << m_FixedImageMax << "]" << std::endl;
    std::cout << "  Moving image range: [" << m_MovingImageMin << ", " << m_MovingImageMax << "]" << std::endl;
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
    
    // 收集所有有效点
    std::vector<ImageType::IndexType> allIndices;
    IteratorType it(m_FixedImage, region);
    
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        allIndices.push_back(it.GetIndex());
    }

    // 随机打乱并选择采样点
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(allIndices.begin(), allIndices.end(), gen);

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

void MattesMutualInformation::ComputeJointPDF(const ParametersType& parameters)
{
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
        // 应用变换
        ImageType::PointType transformedPoint;
        ApplyTransform(sample.point, parameters, transformedPoint);

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

double MattesMutualInformation::GetValue(const ParametersType& parameters)
{
    ComputeJointPDF(parameters);
    double mi = ComputeMutualInformation();
    return -mi; // 返回负值,因为我们要最小化
}

void MattesMutualInformation::GetDerivative(const ParametersType& parameters, 
                                           ParametersType& derivative)
{
    // 使用有限差分法计算梯度
    derivative.resize(parameters.size());
    const double delta = 1e-4;

    double baseValue = GetValue(parameters);

    for (size_t i = 0; i < parameters.size(); ++i)
    {
        ParametersType perturbedParams = parameters;
        perturbedParams[i] += delta;
        double perturbedValue = GetValue(perturbedParams);
        derivative[i] = (perturbedValue - baseValue) / delta;
    }
}

void MattesMutualInformation::GetValueAndDerivative(const ParametersType& parameters,
                                                   double& value,
                                                   ParametersType& derivative)
{
    value = GetValue(parameters);
    GetDerivative(parameters, derivative);
}

void MattesMutualInformation::ApplyTransform(const ImageType::PointType& inputPoint,
                                            const ParametersType& parameters,
                                            ImageType::PointType& outputPoint)
{
    // 简化的仿射变换: 6个参数 (3个平移 + 3个旋转)
    // 实际应用中可以扩展为完整的12参数仿射变换
    
    if (parameters.size() < 6)
    {
        outputPoint = inputPoint;
        return;
    }

    // 提取参数
    double tx = parameters[0];
    double ty = parameters[1];
    double tz = parameters[2];
    double rx = parameters[3]; // 绕X轴旋转(弧度)
    double ry = parameters[4]; // 绕Y轴旋转
    double rz = parameters[5]; // 绕Z轴旋转

    // 计算旋转矩阵 (ZYX欧拉角)
    double cx = std::cos(rx), sx = std::sin(rx);
    double cy = std::cos(ry), sy = std::sin(ry);
    double cz = std::cos(rz), sz = std::sin(rz);

    // 组合旋转矩阵
    double r00 = cy * cz;
    double r01 = -cy * sz;
    double r02 = sy;
    double r10 = sx * sy * cz + cx * sz;
    double r11 = -sx * sy * sz + cx * cz;
    double r12 = -sx * cy;
    double r20 = -cx * sy * cz + sx * sz;
    double r21 = cx * sy * sz + sx * cz;
    double r22 = cx * cy;

    // 应用旋转和平移
    double x = inputPoint[0];
    double y = inputPoint[1];
    double z = inputPoint[2];

    outputPoint[0] = r00 * x + r01 * y + r02 * z + tx;
    outputPoint[1] = r10 * x + r11 * y + r12 * z + ty;
    outputPoint[2] = r20 * x + r21 * y + r22 * z + tz;
}

int MattesMutualInformation::ComputeFixedImageBin(double value) const
{
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
    if (value <= m_MovingImageMin)
        return 0;
    if (value >= m_MovingImageMax)
        return m_NumberOfHistogramBins - 1;

    double range = m_MovingImageMax - m_MovingImageMin;
    int bin = static_cast<int>((value - m_MovingImageMin) / range * m_NumberOfHistogramBins);
    
    return std::min(bin, static_cast<int>(m_NumberOfHistogramBins - 1));
}

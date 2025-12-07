#include "ImageRegistration.h"
#include "MattesMutualInformation.h"
#include "RegularStepGradientDescentOptimizer.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <itkImageFileReader.h>

ImageRegistration::ImageRegistration()
    : m_NumberOfHistogramBins(50)
    , m_NumberOfSpatialSamples(100000)
    , m_LearningRate(0.5)
    , m_MinimumStepLength(0.0001)
    , m_NumberOfIterations(300)
    , m_RelaxationFactor(0.8)
    , m_GradientMagnitudeTolerance(1e-4)
    , m_NumberOfLevels(3)
    , m_RandomSeed(121212)
    , m_FinalMetricValue(0.0)
    , m_ElapsedTime(0.0)
{
    // 默认多分辨率设置
    m_ShrinkFactors = {4, 2, 1};
    m_SmoothingSigmas = {2.0, 1.0, 0.0};
    
    m_Metric = std::make_unique<MattesMutualInformation>();
    m_Optimizer = std::make_unique<RegularStepGradientDescentOptimizer>();
    m_Transform = TransformType::New();
}

ImageRegistration::~ImageRegistration()
{
}

void ImageRegistration::SetFixedImagePath(const std::string& path)
{
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName(path);
    try
    {
        reader->Update();
        m_FixedImage = reader->GetOutput();
        std::cout << "  Fixed image loaded: " << path << std::endl;
    }
    catch (const itk::ExceptionObject& e)
    {
        std::cerr << "Error loading fixed image: " << e.what() << std::endl;
        throw;
    }
}

void ImageRegistration::SetMovingImagePath(const std::string& path)
{
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName(path);
    try
    {
        reader->Update();
        m_MovingImage = reader->GetOutput();
        std::cout << "  Moving image loaded: " << path << std::endl;
    }
    catch (const itk::ExceptionObject& e)
    {
        std::cerr << "Error loading moving image: " << e.what() << std::endl;
        throw;
    }
}

void ImageRegistration::SetFixedImage(ImageType::Pointer fixedImage)
{
    m_FixedImage = fixedImage;
}

void ImageRegistration::SetMovingImage(ImageType::Pointer movingImage)
{
    m_MovingImage = movingImage;
}

void ImageRegistration::ComputeGeometricCenter(ImageType::Pointer image, ImageType::PointType& center)
{
    auto region = image->GetLargestPossibleRegion();
    auto size = region.GetSize();
    auto spacing = image->GetSpacing();
    auto origin = image->GetOrigin();
    auto direction = image->GetDirection();

    // 计算图像中心索引
    ImageType::IndexType centerIndex;
    for (unsigned int i = 0; i < 3; ++i)
    {
        centerIndex[i] = size[i] / 2;
    }

    // 转换为物理坐标
    image->TransformIndexToPhysicalPoint(centerIndex, center);
}

void ImageRegistration::InitializeTransform()
{
    // 使用几何中心初始化变换
    ImageType::PointType fixedCenter, movingCenter;
    ComputeGeometricCenter(m_FixedImage, fixedCenter);
    ComputeGeometricCenter(m_MovingImage, movingCenter);

    // 设置旋转中心为固定图像的几何中心
    m_Transform->SetCenter(fixedCenter);
    
    // 计算初始平移(使两个图像中心对齐)
    TransformType::OutputVectorType initialTranslation;
    initialTranslation[0] = movingCenter[0] - fixedCenter[0];
    initialTranslation[1] = movingCenter[1] - fixedCenter[1];
    initialTranslation[2] = movingCenter[2] - fixedCenter[2];
    m_Transform->SetTranslation(initialTranslation);
    
    // 初始旋转为0
    m_Transform->SetRotation(0.0, 0.0, 0.0);

    std::cout << "Initial transform center: [" 
              << fixedCenter[0] << ", " << fixedCenter[1] << ", " << fixedCenter[2] << "]" << std::endl;
}

ImageRegistration::ImageType::Pointer ImageRegistration::ShrinkImage(ImageType::Pointer image, unsigned int factor)
{
    if (factor <= 1)
    {
        return image;
    }

    using ShrinkFilterType = itk::ShrinkImageFilter<ImageType, ImageType>;
    auto shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetInput(image);
    shrinkFilter->SetShrinkFactors(factor);
    shrinkFilter->Update();
    return shrinkFilter->GetOutput();
}

ImageRegistration::ImageType::Pointer ImageRegistration::SmoothImage(ImageType::Pointer image, double sigma)
{
    if (sigma <= 0.0)
    {
        return image;
    }

    using SmoothingFilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
    auto smoothFilter = SmoothingFilterType::New();
    smoothFilter->SetInput(image);
    smoothFilter->SetSigma(sigma);
    smoothFilter->Update();
    return smoothFilter->GetOutput();
}

void ImageRegistration::RunSingleLevel(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level)
{
    // 配置度量
    m_Metric->SetFixedImage(fixedImage);
    m_Metric->SetMovingImage(movingImage);
    m_Metric->SetNumberOfHistogramBins(m_NumberOfHistogramBins);
    m_Metric->SetNumberOfSpatialSamples(m_NumberOfSpatialSamples);
    m_Metric->SetRandomSeed(m_RandomSeed);
    m_Metric->SetTransform(m_Transform);
    m_Metric->SetUseStratifiedSampling(true);  // 使用分层采样
    m_Metric->Initialize();

    // 配置优化器
    m_Optimizer->SetLearningRate(m_LearningRate);
    m_Optimizer->SetMinimumStepLength(m_MinimumStepLength);
    m_Optimizer->SetNumberOfIterations(m_NumberOfIterations);
    m_Optimizer->SetRelaxationFactor(m_RelaxationFactor);
    m_Optimizer->SetGradientMagnitudeTolerance(m_GradientMagnitudeTolerance);
    m_Optimizer->SetReturnBestParametersAndValue(true);
    m_Optimizer->SetNumberOfParameters(6);
    
    // 使用自动估算的参数尺度
    std::vector<double> scales = EstimateParameterScales();
    m_Optimizer->SetScales(scales);

    // 设置代价函数
    m_Optimizer->SetCostFunction([this]() -> double {
        return m_Metric->GetValue();
    });

    // 设置梯度函数
    m_Optimizer->SetGradientFunction([this](std::vector<double>& gradient) {
        m_Metric->GetDerivative(gradient);
    });
    
    // 设置参数获取函数
    m_Optimizer->SetGetParametersFunction([this]() -> std::vector<double> {
        auto params = m_Transform->GetParameters();
        std::vector<double> result(params.Size());
        for (unsigned int i = 0; i < params.Size(); ++i)
        {
            result[i] = params[i];
        }
        return result;
    });
    
    // 设置参数设置函数
    m_Optimizer->SetSetParametersFunction([this](const std::vector<double>& params) {
        TransformType::ParametersType itkParams(6);
        for (unsigned int i = 0; i < 6 && i < params.size(); ++i)
        {
            itkParams[i] = params[i];
        }
        m_Transform->SetParameters(itkParams);
    });

    // 设置参数更新函数 (不再使用，保留兼容性)
    m_Optimizer->SetUpdateParametersFunction([this](const std::vector<double>& update) {
        auto params = m_Transform->GetParameters();
        for (unsigned int i = 0; i < 6 && i < params.Size(); ++i)
        {
            params[i] += update[i];
        }
        m_Transform->SetParameters(params);
    });

    // 设置观察者输出
    if (m_IterationObserver)
    {
        m_Optimizer->SetObserver([this](unsigned int iter, double value, double stepLength) {
            m_IterationObserver(iter, value, stepLength);
        });
    }
    else
    {
        m_Optimizer->SetObserver([](unsigned int iter, double value, double stepLength) {
            std::cout << "  Iter: " << std::setw(4) << iter 
                      << "  Metric: " << std::setw(12) << std::fixed << std::setprecision(6) << value
                      << "  LearningRate: " << std::setw(10) << std::scientific << std::setprecision(4) << stepLength
                      << std::endl;
        });
    }

    // 执行优化
    m_Optimizer->StartOptimization();

    m_FinalMetricValue = m_Optimizer->GetBestValue();
}

void ImageRegistration::Update()
{
    if (!m_FixedImage || !m_MovingImage)
    {
        throw std::runtime_error("Fixed or moving image not set");
    }

    auto startTime = std::chrono::high_resolution_clock::now();

    // 初始化变换
    InitializeTransform();

    // 打印多分辨率策略
    std::cout << "\nMulti-Resolution Strategy: " << m_NumberOfLevels << " levels" << std::endl;
    for (unsigned int i = 0; i < m_NumberOfLevels; ++i)
    {
        std::cout << "  Level " << i << ": Shrink " << m_ShrinkFactors[i] 
                  << "x, Smooth " << m_SmoothingSigmas[i] << " mm";
        if (i == 0) std::cout << " (coarse)";
        else if (i == m_NumberOfLevels - 1) std::cout << " (fine)";
        else std::cout << " (medium)";
        std::cout << std::endl;
    }

    // 计算采样信息
    auto fixedRegion = m_FixedImage->GetLargestPossibleRegion();
    auto fixedSize = fixedRegion.GetSize();
    unsigned long totalVoxels = fixedSize[0] * fixedSize[1] * fixedSize[2];
    double samplingPercentage = (double)m_NumberOfSpatialSamples / (double)totalVoxels * 100.0;
    std::cout << "Using " << m_NumberOfSpatialSamples << " samples (" 
              << std::fixed << std::setprecision(1) << samplingPercentage 
              << "% of total voxels)\n";

    // 多分辨率金字塔
    for (unsigned int level = 0; level < m_NumberOfLevels; ++level)
    {
        // 调用级别观察者
        unsigned int shrinkFactor = (level < m_ShrinkFactors.size()) ? m_ShrinkFactors[level] : 1;
        double smoothingSigma = (level < m_SmoothingSigmas.size()) ? m_SmoothingSigmas[level] : 0.0;
        
        if (m_LevelObserver)
        {
            m_LevelObserver(level, shrinkFactor, smoothingSigma);
        }
        else
        {
            std::cout << "\n=== Multi-Resolution Level " << level << " of " << (m_NumberOfLevels - 1) << " ===" << std::endl;
        }

        // 下采样和平滑
        ImageType::Pointer fixedPyramid = ShrinkImage(m_FixedImage, shrinkFactor);
        fixedPyramid = SmoothImage(fixedPyramid, smoothingSigma);

        ImageType::Pointer movingPyramid = ShrinkImage(m_MovingImage, shrinkFactor);
        movingPyramid = SmoothImage(movingPyramid, smoothingSigma);

        // 在当前分辨率运行配准
        RunSingleLevel(fixedPyramid, movingPyramid, level);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    m_ElapsedTime = std::chrono::duration<double>(endTime - startTime).count();

    std::cout << "Final metric value: " << std::scientific << std::setprecision(4) 
              << m_FinalMetricValue << std::endl;
}

// ============================================================================
// 计算图像的物理半径 (用于参数尺度估算)
// ============================================================================

double ImageRegistration::ComputePhysicalRadius(ImageType::Pointer image)
{
    // 计算图像对角线长度的一半作为"物理半径"
    auto region = image->GetLargestPossibleRegion();
    auto size = region.GetSize();
    auto spacing = image->GetSpacing();
    
    double diagonalSquared = 0.0;
    for (unsigned int dim = 0; dim < 3; ++dim)
    {
        double physicalLength = size[dim] * spacing[dim];
        diagonalSquared += physicalLength * physicalLength;
    }
    
    return std::sqrt(diagonalSquared) / 2.0;
}

// ============================================================================
// 自动估算参数尺度 (类似ITK的RegistrationParameterScalesFromPhysicalShift)
// ============================================================================

std::vector<double> ImageRegistration::EstimateParameterScales()
{
    // ITK的RegistrationParameterScalesFromPhysicalShift原理:
    // 对于每个参数,估算单位参数变化导致的物理空间平均位移
    // scales[i] = 1.0 / average_physical_shift[i]
    // 这样让所有参数在优化时有相似的"步幅"
    
    std::vector<double> scales(6);
    
    // 计算图像的物理半径 (用于估算旋转的影响)
    double physicalRadius = ComputePhysicalRadius(m_FixedImage);
    
    // 旋转参数的尺度估算:
    // 旋转1弧度在半径R处产生的位移约为R
    // 所以 scale_rotation = 1.0 / R
    // 这使得旋转参数的"有效步长"与平移参数可比
    double rotationScale = 1.0 / physicalRadius;
    
    // 平移参数的尺度:
    // 平移1mm产生的位移就是1mm
    // scale_translation = 1.0 / 1.0 = 1.0
    double translationScale = 1.0;
    
    scales[0] = rotationScale;  // rotation X
    scales[1] = rotationScale;  // rotation Y
    scales[2] = rotationScale;  // rotation Z
    scales[3] = translationScale;  // translation X
    scales[4] = translationScale;  // translation Y
    scales[5] = translationScale;  // translation Z
    
    std::cout << "Parameter Scales (auto-estimated):" << std::endl;
    std::cout << "  Physical radius: " << std::fixed << std::setprecision(2) 
              << physicalRadius << " mm" << std::endl;
    std::cout << "  Rotation scale: " << std::scientific << std::setprecision(4) 
              << rotationScale << std::endl;
    std::cout << "  Translation scale: " << std::fixed << std::setprecision(4) 
              << translationScale << std::endl;
    
    return scales;
}

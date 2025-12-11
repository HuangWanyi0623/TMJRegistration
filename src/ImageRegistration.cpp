#include "ImageRegistration.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <itkImageFileReader.h>
#include <itkTransformFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

// ============================================================================
// 构造函数和析构函数
// ============================================================================

ImageRegistration::ImageRegistration()
    : m_TransformType(ConfigManager::TransformType::Rigid)
    , m_UseInitialTransform(false)
    , m_NumberOfHistogramBins(50)
    , m_NumberOfSpatialSamples(100000)
    , m_LearningRate(0.5)
    , m_MinimumStepLength(0.0001)
    , m_NumberOfIterations({300})  // 默认单层300次迭代
    , m_RelaxationFactor(0.8)
    , m_GradientMagnitudeTolerance(1e-4)
    , m_NumberOfLevels(3)
    , m_RandomSeed(121212)
    , m_UseStratifiedSampling(true)
    , m_SamplingPercentage(0.10)
    , m_FinalMetricValue(0.0)
    , m_ElapsedTime(0.0)
    , m_Verbose(false)
    , m_MaskVoxelCount(0)  // 初始化掩膜体素数为0
{
    // 默认多分辨率设置
    m_ShrinkFactors = {4, 2, 1};
    m_SmoothingSigmas = {2.0, 1.0, 0.0};
    
    m_Metric = std::make_unique<MattesMutualInformation>();
    m_Optimizer = std::make_unique<RegularStepGradientDescentOptimizer>();
    
    // 初始化变换
    m_RigidTransform = RigidTransformType::New();
    m_AffineTransform = AffineTransformType::New();
    m_InitialTransform = CompositeTransformType::New();
}

ImageRegistration::~ImageRegistration()
{
}

// ============================================================================
// 图像加载
// ============================================================================

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

// ============================================================================
// 掩膜加载 (用于局部配准)
// ============================================================================

bool ImageRegistration::LoadFixedMask(const std::string& maskFilePath)
{
    try
    {
        // 读取掩膜图像 (通常是 unsigned char 类型的 LabelMap)
        using MaskReaderType = itk::ImageFileReader<MaskImageType>;
        auto maskReader = MaskReaderType::New();
        maskReader->SetFileName(maskFilePath);
        maskReader->Update();
        
        MaskImageType::Pointer maskImage = maskReader->GetOutput();
        
        // 将图像包装为 ImageMaskSpatialObject
        m_FixedImageMask = MaskSpatialObjectType::New();
        m_FixedImageMask->SetImage(maskImage);
        m_FixedImageMask->Update();
        
        // 统计掩膜内的体素数量
        using IteratorType = itk::ImageRegionConstIterator<MaskImageType>;
        IteratorType it(maskImage, maskImage->GetLargestPossibleRegion());
        unsigned long maskVoxels = 0;
        unsigned long totalVoxels = 0;
        for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            ++totalVoxels;
            if (it.Get() > 0)
            {
                ++maskVoxels;
            }
        }
        
        double maskPercentage = 100.0 * static_cast<double>(maskVoxels) / static_cast<double>(totalVoxels);
        
        // 保存掩膜体素数供后续使用
        m_MaskVoxelCount = maskVoxels;
        
        std::cout << "[Fixed Mask] Loaded: " << maskFilePath << std::endl;
        std::cout << "  Mask coverage: " << maskVoxels << " / " << totalVoxels 
                  << " voxels (" << std::fixed << std::setprecision(1) << maskPercentage << "%)" << std::endl;
        
        return true;
    }
    catch (const itk::ExceptionObject& e)
    {
        std::cerr << "[Error] Failed to load mask: " << e.what() << std::endl;
        m_FixedImageMask = nullptr;
        return false;
    }
}

// ============================================================================
// 变换类型设置
// ============================================================================

void ImageRegistration::SetTransformType(ConfigManager::TransformType type)
{
    m_TransformType = type;
}

unsigned int ImageRegistration::GetNumberOfParameters() const
{
    return (m_TransformType == ConfigManager::TransformType::Rigid) ? 6 : 12;
}

// ============================================================================
// 初始变换加载
// ============================================================================

bool ImageRegistration::LoadInitialTransform(const std::string& h5FilePath)
{
    try
    {
        auto reader = itk::TransformFileReader::New();
        reader->SetFileName(h5FilePath);
        reader->Update();
        
        auto transformList = reader->GetTransformList();
        if (transformList->empty())
        {
            std::cerr << "[Warning] No transforms found in: " << h5FilePath << std::endl;
            return false;
        }
        
        // 清空现有的初始变换
        m_InitialTransform = CompositeTransformType::New();
        
        // 添加所有读取的变换到复合变换
        // 注意：ITK的CompositeTransform是后进先出（LIFO），最后添加的变换最先应用
        for (auto it = transformList->begin(); it != transformList->end(); ++it)
        {
            // 尝试转换为3D变换
            auto transform3D = dynamic_cast<itk::Transform<double, 3, 3>*>((*it).GetPointer());
            if (transform3D)
            {
                m_InitialTransform->AddTransform(transform3D);
            }
        }
        
        m_UseInitialTransform = true;
        std::cout << "[Initial Transform] Loaded from: " << h5FilePath << std::endl;
        std::cout << "  Number of transforms: " << m_InitialTransform->GetNumberOfTransforms() << std::endl;
        
        // 调试：测试初始变换是否使图像对齐
        // 取Fixed图像中心点，看变换后是否接近Moving图像中心
        ImageType::PointType fixedCenter, movingCenter;
        ComputeGeometricCenter(m_FixedImage, fixedCenter);
        ComputeGeometricCenter(m_MovingImage, movingCenter);
        
        ImageType::PointType transformedFixedCenter = m_InitialTransform->TransformPoint(fixedCenter);
        
        double distance = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            double diff = transformedFixedCenter[i] - movingCenter[i];
            distance += diff * diff;
        }
        distance = std::sqrt(distance);
        
        std::cout << "[Initial Transform Verification]" << std::endl;
        std::cout << "  Fixed center: [" << fixedCenter[0] << ", " << fixedCenter[1] << ", " << fixedCenter[2] << "]" << std::endl;
        std::cout << "  Moving center: [" << movingCenter[0] << ", " << movingCenter[1] << ", " << movingCenter[2] << "]" << std::endl;
        std::cout << "  Transformed fixed center: [" << transformedFixedCenter[0] << ", " << transformedFixedCenter[1] << ", " << transformedFixedCenter[2] << "]" << std::endl;
        std::cout << "  Distance to moving center: " << std::fixed << std::setprecision(2) << distance << " mm" << std::endl;
        
        if (distance > 50.0)  // 如果距离超过50mm，可能变换方向错了
        {
            std::cout << "  [WARNING] Initial transform may be in wrong direction! Distance > 50mm" << std::endl;
            std::cout << "  [HINT] If registration fails, try using inverse transform or check transform direction." << std::endl;
        }
        else
        {
            std::cout << "  [OK] Initial transform appears correct (distance < 50mm)" << std::endl;
        }
        
        return true;
    }
    catch (const itk::ExceptionObject& e)
    {
        std::cerr << "[Error] Failed to load initial transform: " << e.what() << std::endl;
        return false;
    }
}

// ============================================================================
// 从配置加载
// ============================================================================

void ImageRegistration::LoadFromConfig(const ConfigManager::RegistrationConfig& config)
{
    m_TransformType = config.transformType;
    m_NumberOfHistogramBins = config.numberOfHistogramBins;
    m_NumberOfSpatialSamples = config.numberOfSpatialSamples;
    m_LearningRate = config.learningRate;
    m_MinimumStepLength = config.minimumStepLength;
    m_NumberOfIterations = config.numberOfIterations;  // 现在是vector
    m_RelaxationFactor = config.relaxationFactor;
    m_GradientMagnitudeTolerance = config.gradientMagnitudeTolerance;
    m_NumberOfLevels = config.numberOfLevels;
    m_ShrinkFactors = config.shrinkFactors;
    m_SmoothingSigmas = config.smoothingSigmas;
    m_UseStratifiedSampling = config.useStratifiedSampling;
    m_RandomSeed = config.randomSeed;
    m_SamplingPercentage = config.samplingPercentage;
}

// ============================================================================
// 获取复合变换结果
// ============================================================================

ImageRegistration::CompositeTransformType::Pointer ImageRegistration::GetFinalCompositeTransform() const
{
    auto composite = CompositeTransformType::New();
    
    // 关键修复：如果使用了初始变换并且已经将其参数初始化到当前变换中，
    // 那么当前变换已经包含了所有信息，不应该再添加初始变换！
    // 这避免了粗配准被应用两次的问题。
    
    // 只添加优化后的最终变换
    if (m_TransformType == ConfigManager::TransformType::Rigid)
    {
        composite->AddTransform(m_RigidTransform);
    }
    else
    {
        composite->AddTransform(m_AffineTransform);
    }
    
    return composite;
}

// ============================================================================
// 几何中心计算
// ============================================================================

void ImageRegistration::ComputeGeometricCenter(ImageType::Pointer image, ImageType::PointType& center)
{
    auto region = image->GetLargestPossibleRegion();
    auto size = region.GetSize();

    ImageType::IndexType centerIndex;
    for (unsigned int i = 0; i < 3; ++i)
    {
        centerIndex[i] = size[i] / 2;
    }

    image->TransformIndexToPhysicalPoint(centerIndex, center);
}

// ============================================================================
// 变换初始化
// ============================================================================

void ImageRegistration::InitializeTransform()
{
    if (m_TransformType == ConfigManager::TransformType::Rigid)
    {
        InitializeRigidTransform();
    }
    else
    {
        InitializeAffineTransform();
    }
}

void ImageRegistration::InitializeRigidTransform()
{
    ImageType::PointType fixedCenter, movingCenter;
    ComputeGeometricCenter(m_FixedImage, fixedCenter);
    ComputeGeometricCenter(m_MovingImage, movingCenter);

    // 如果有初始变换，从中提取参数初始化当前变换
    // 关键：必须保留初始变换的旋转中心！否则旋转角度相同但效果完全不同
    if (m_UseInitialTransform && m_InitialTransform->GetNumberOfTransforms() > 0)
    {
        std::cout << "[Initial Transform] Initializing current transform from loaded transform..." << std::endl;
        
        // 尝试获取初始变换的参数
        auto firstTransform = m_InitialTransform->GetNthTransform(0);
        
        if (auto rigidInit = dynamic_cast<RigidTransformType*>(firstTransform.GetPointer()))
        {
            // 如果初始变换也是刚体，直接复制所有参数（包括中心点！）
            m_RigidTransform->SetCenter(rigidInit->GetCenter());  // <- 关键：保留原始中心
            m_RigidTransform->SetParameters(rigidInit->GetParameters());
            
            std::cout << "  Initialized from Rigid transform parameters" << std::endl;
            auto params = rigidInit->GetParameters();
            auto center = rigidInit->GetCenter();
            std::cout << "    Rotation center: [" << center[0] << ", " << center[1] << ", " << center[2] << "]" << std::endl;
            std::cout << "    Rotation (rad): [" << params[0] << ", " << params[1] << ", " << params[2] << "]" << std::endl;
            std::cout << "    Translation (mm): [" << params[3] << ", " << params[4] << ", " << params[5] << "]" << std::endl;
        }
        else if (auto affineInit = dynamic_cast<AffineTransformType*>(firstTransform.GetPointer()))
        {
            // 如果初始变换是仿射，保留其中心点并提取旋转和平移部分
            m_RigidTransform->SetCenter(affineInit->GetCenter());  // <- 关键：保留原始中心
            
            auto matrix = affineInit->GetMatrix();
            auto translation = affineInit->GetTranslation();
            
            // 从仿射矩阵中提取欧拉角（简化方法，假设矩阵接近正交）
            double angleX = std::atan2(matrix(2, 1), matrix(2, 2));
            double angleY = std::atan2(-matrix(2, 0), std::sqrt(matrix(2, 1) * matrix(2, 1) + matrix(2, 2) * matrix(2, 2)));
            double angleZ = std::atan2(matrix(1, 0), matrix(0, 0));
            
            m_RigidTransform->SetRotation(angleX, angleY, angleZ);
            m_RigidTransform->SetTranslation(translation);
            
            std::cout << "  Initialized from Affine transform (extracted rigid component)" << std::endl;
            auto center = affineInit->GetCenter();
            std::cout << "    Rotation center: [" << center[0] << ", " << center[1] << ", " << center[2] << "]" << std::endl;
            std::cout << "    Rotation (rad): [" << angleX << ", " << angleY << ", " << angleZ << "]" << std::endl;
            std::cout << "    Translation (mm): [" << translation[0] << ", " << translation[1] << ", " << translation[2] << "]" << std::endl;
        }
        else
        {
            std::cout << "  [Warning] Initial transform type not recognized, using default initialization" << std::endl;
            // 回退到默认初始化
            m_RigidTransform->SetCenter(fixedCenter);
            RigidTransformType::OutputVectorType initialTranslation;
            initialTranslation[0] = movingCenter[0] - fixedCenter[0];
            initialTranslation[1] = movingCenter[1] - fixedCenter[1];
            initialTranslation[2] = movingCenter[2] - fixedCenter[2];
            m_RigidTransform->SetTranslation(initialTranslation);
            m_RigidTransform->SetRotation(0.0, 0.0, 0.0);
        }
    }
    else
    {
        // 没有初始变换，使用默认初始化（图像中心对齐）
        m_RigidTransform->SetCenter(fixedCenter);
        RigidTransformType::OutputVectorType initialTranslation;
        initialTranslation[0] = movingCenter[0] - fixedCenter[0];
        initialTranslation[1] = movingCenter[1] - fixedCenter[1];
        initialTranslation[2] = movingCenter[2] - fixedCenter[2];
        m_RigidTransform->SetTranslation(initialTranslation);
        m_RigidTransform->SetRotation(0.0, 0.0, 0.0);
        
        std::cout << "Initial transform center: [" 
                  << std::fixed << std::setprecision(2)
                  << fixedCenter[0] << ", " << fixedCenter[1] << ", " << fixedCenter[2] << "]" << std::endl;
    }
}

void ImageRegistration::InitializeAffineTransform()
{
    ImageType::PointType fixedCenter, movingCenter;
    ComputeGeometricCenter(m_FixedImage, fixedCenter);
    ComputeGeometricCenter(m_MovingImage, movingCenter);

    // 如果有初始变换(例如刚体配准结果),从中初始化仿射变换
    if (m_UseInitialTransform && m_InitialTransform->GetNumberOfTransforms() > 0)
    {
        std::cout << "[Initial Transform] Initializing Affine from loaded transform..." << std::endl;
        
        auto firstTransform = m_InitialTransform->GetNthTransform(0);
        
        if (auto rigidInit = dynamic_cast<RigidTransformType*>(firstTransform.GetPointer()))
        {
            // 从刚体变换初始化:保留中心、矩阵和平移
            m_AffineTransform->SetCenter(rigidInit->GetCenter());
            m_AffineTransform->SetMatrix(rigidInit->GetMatrix());
            m_AffineTransform->SetTranslation(rigidInit->GetTranslation());
            
            std::cout << "  Initialized from Rigid transform" << std::endl;
            auto center = rigidInit->GetCenter();
            auto params = rigidInit->GetParameters();
            std::cout << "    Center: [" << center[0] << ", " << center[1] << ", " << center[2] << "]" << std::endl;
            std::cout << "    Rotation (rad): [" << params[0] << ", " << params[1] << ", " << params[2] << "]" << std::endl;
            std::cout << "    Translation (mm): [" << params[3] << ", " << params[4] << ", " << params[5] << "]" << std::endl;
        }
        else if (auto affineInit = dynamic_cast<AffineTransformType*>(firstTransform.GetPointer()))
        {
            // 从仿射变换初始化:直接复制
            m_AffineTransform->SetCenter(affineInit->GetCenter());
            m_AffineTransform->SetMatrix(affineInit->GetMatrix());
            m_AffineTransform->SetTranslation(affineInit->GetTranslation());
            
            std::cout << "  Initialized from Affine transform" << std::endl;
            auto center = affineInit->GetCenter();
            std::cout << "    Center: [" << center[0] << ", " << center[1] << ", " << center[2] << "]" << std::endl;
        }
        else
        {
            std::cout << "  [Warning] Initial transform type not recognized, using default initialization" << std::endl;
            // 回退到默认初始化
            m_AffineTransform->SetCenter(fixedCenter);
            AffineTransformType::MatrixType matrix;
            matrix.SetIdentity();
            m_AffineTransform->SetMatrix(matrix);
            
            AffineTransformType::OutputVectorType initialTranslation;
            initialTranslation[0] = movingCenter[0] - fixedCenter[0];
            initialTranslation[1] = movingCenter[1] - fixedCenter[1];
            initialTranslation[2] = movingCenter[2] - fixedCenter[2];
            m_AffineTransform->SetTranslation(initialTranslation);
        }
    }
    else
    {
        // 没有初始变换,使用默认初始化
        m_AffineTransform->SetCenter(fixedCenter);
        
        AffineTransformType::MatrixType matrix;
        matrix.SetIdentity();
        m_AffineTransform->SetMatrix(matrix);
        
        AffineTransformType::OutputVectorType initialTranslation;
        initialTranslation[0] = movingCenter[0] - fixedCenter[0];
        initialTranslation[1] = movingCenter[1] - fixedCenter[1];
        initialTranslation[2] = movingCenter[2] - fixedCenter[2];
        m_AffineTransform->SetTranslation(initialTranslation);

        std::cout << "Initial affine transform center: [" 
                  << fixedCenter[0] << ", " << fixedCenter[1] << ", " << fixedCenter[2] << "]" << std::endl;
    }
}

// ============================================================================
// 图像预处理
// ============================================================================

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

ImageRegistration::ImageType::Pointer ImageRegistration::WinsorizeImage(ImageType::Pointer image, double lowerQuantile, double upperQuantile)
{
    // ANTs风格的Winsorizing: 将强度截断在指定分位数之间
    // 这能极大提高MI对软组织的敏感度,忽略骨骼高亮或背景噪声
    
    using IteratorType = itk::ImageRegionConstIterator<ImageType>;
    using OutputIteratorType = itk::ImageRegionIterator<ImageType>;
    
    // 1. 收集所有像素值
    std::vector<float> values;
    IteratorType it(image, image->GetLargestPossibleRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        values.push_back(it.Get());
    }
    
    // 2. 排序并计算分位数
    std::sort(values.begin(), values.end());
    size_t lowerIndex = static_cast<size_t>(lowerQuantile * values.size());
    size_t upperIndex = static_cast<size_t>(upperQuantile * values.size());
    
    if (lowerIndex >= values.size()) lowerIndex = 0;
    if (upperIndex >= values.size()) upperIndex = values.size() - 1;
    
    float lowerThreshold = values[lowerIndex];
    float upperThreshold = values[upperIndex];
    
    // 3. 创建输出图像并截断
    auto output = ImageType::New();
    output->SetRegions(image->GetLargestPossibleRegion());
    output->SetSpacing(image->GetSpacing());
    output->SetOrigin(image->GetOrigin());
    output->SetDirection(image->GetDirection());
    output->Allocate();
    
    IteratorType inputIt(image, image->GetLargestPossibleRegion());
    OutputIteratorType outputIt(output, output->GetLargestPossibleRegion());
    
    for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
    {
        float value = inputIt.Get();
        if (value < lowerThreshold)
            value = lowerThreshold;
        else if (value > upperThreshold)
            value = upperThreshold;
        outputIt.Set(value);
    }
    
    return output;
}

// ============================================================================
// 参数尺度估算
// ============================================================================

double ImageRegistration::ComputePhysicalRadius(ImageType::Pointer image)
{
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

std::vector<double> ImageRegistration::EstimateParameterScales()
{
    if (m_TransformType == ConfigManager::TransformType::Rigid)
    {
        return EstimateRigidParameterScales();
    }
    else
    {
        return EstimateAffineParameterScales();
    }
}

std::vector<double> ImageRegistration::EstimateRigidParameterScales()
{
    std::vector<double> scales(6);
    double physicalRadius = ComputePhysicalRadius(m_FixedImage);
    
    // rotationScale should be in physical units: larger radius -> larger rotation effect in mm
    double rotationScale = physicalRadius;
    double translationScale = 1.0;
    
    scales[0] = rotationScale;
    scales[1] = rotationScale;
    scales[2] = rotationScale;
    scales[3] = translationScale;
    scales[4] = translationScale;
    scales[5] = translationScale;
    
    std::cout << "Parameter Scales (auto-estimated):" << std::endl;
    std::cout << "  Physical radius: " << std::fixed << std::setprecision(2) 
              << physicalRadius << " mm" << std::endl;
    std::cout << "  Rotation scale: " << std::fixed << std::setprecision(2) 
              << rotationScale << " (physical radius)" << std::endl;
    std::cout << "  Translation scale: " << std::fixed << std::setprecision(4) 
              << translationScale << std::endl;
    
    return scales;
}

std::vector<double> ImageRegistration::EstimateAffineParameterScales()
{
    // 仿射变换有12个参数: 9个矩阵元素 + 3个平移
    // AffineTransform参数顺序: [M00, M01, M02, M10, M11, M12, M20, M21, M22, Tx, Ty, Tz]
    std::vector<double> scales(12);
    double physicalRadius = ComputePhysicalRadius(m_FixedImage);
    
    // 矩阵元素对位移的影响约为 radius
    double matrixScale = physicalRadius;
    double translationScale = 1.0;
    
    for (int i = 0; i < 9; ++i)
    {
        scales[i] = matrixScale;
    }
    for (int i = 9; i < 12; ++i)
    {
        scales[i] = translationScale;
    }
    
    std::cout << "Parameter Scales (Affine, auto-estimated):" << std::endl;
    std::cout << "  Physical radius: " << std::fixed << std::setprecision(2) 
              << physicalRadius << " mm" << std::endl;
    std::cout << "  Matrix scale: " << std::fixed << std::setprecision(2) 
              << matrixScale << " (physical radius)" << std::endl;
    std::cout << "  Translation scale: " << std::fixed << std::setprecision(4) 
              << translationScale << std::endl;
    
    return scales;
}

// ============================================================================
// 单层配准 - 主分发函数
// ============================================================================

void ImageRegistration::RunSingleLevel(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level)
{
    if (m_TransformType == ConfigManager::TransformType::Rigid)
    {
        RunSingleLevelRigid(fixedImage, movingImage, level);
    }
    else
    {
        RunSingleLevelAffine(fixedImage, movingImage, level);
    }
}

// ============================================================================
// 刚体变换的雅可比矩阵计算
// ============================================================================

static void ComputeRigidJacobian(
    const itk::Point<double, 3>& point,
    const itk::Euler3DTransform<double>::Pointer& transform,
    std::vector<std::array<double, 3>>& jacobian)
{
    auto params = transform->GetParameters();
    auto center = transform->GetCenter();
    
    double rx = params[0], ry = params[1], rz = params[2];
    double px = point[0] - center[0];
    double py = point[1] - center[1];
    double pz = point[2] - center[2];
    
    double cx = std::cos(rx), sx = std::sin(rx);
    double cy = std::cos(ry), sy = std::sin(ry);
    double cz = std::cos(rz), sz = std::sin(rz);
    
    jacobian.resize(6);
    
    jacobian[0][0] = 0.0;
    jacobian[0][1] = (-sx*sy*cz - cx*sz) * px + (-sx*sy*sz + cx*cz) * py + (-sx*cy) * pz;
    jacobian[0][2] = (cx*sy*cz - sx*sz) * px + (cx*sy*sz + sx*cz) * py + (cx*cy) * pz;
    
    jacobian[1][0] = (-sy*cz) * px + (-sy*sz) * py + (-cy) * pz;
    jacobian[1][1] = (cx*cy*cz) * px + (cx*cy*sz) * py + (-cx*sy) * pz;
    jacobian[1][2] = (sx*cy*cz) * px + (sx*cy*sz) * py + (-sx*sy) * pz;
    
    jacobian[2][0] = (-cy*sz) * px + (cy*cz) * py;
    jacobian[2][1] = (-sx*sy*sz - cx*cz) * px + (sx*sy*cz - cx*sz) * py;
    jacobian[2][2] = (cx*sy*sz - sx*cz) * px + (-cx*sy*cz - sx*sz) * py;
    
    jacobian[3] = {1.0, 0.0, 0.0};
    jacobian[4] = {0.0, 1.0, 0.0};
    jacobian[5] = {0.0, 0.0, 1.0};
}

// ============================================================================
// 仿射变换的雅可比矩阵计算
// ============================================================================

static void ComputeAffineJacobian(
    const itk::Point<double, 3>& point,
    const itk::AffineTransform<double, 3>::Pointer& transform,
    std::vector<std::array<double, 3>>& jacobian)
{
    // AffineTransform: y = M * (x - center) + center + translation
    // 参数顺序: [M00, M01, M02, M10, M11, M12, M20, M21, M22, Tx, Ty, Tz]
    
    auto center = transform->GetCenter();
    double px = point[0] - center[0];
    double py = point[1] - center[1];
    double pz = point[2] - center[2];
    
    jacobian.resize(12);
    
    // dy/dM00 = [px, 0, 0]
    jacobian[0] = {px, 0.0, 0.0};
    // dy/dM01 = [py, 0, 0]
    jacobian[1] = {py, 0.0, 0.0};
    // dy/dM02 = [pz, 0, 0]
    jacobian[2] = {pz, 0.0, 0.0};
    
    // dy/dM10 = [0, px, 0]
    jacobian[3] = {0.0, px, 0.0};
    // dy/dM11 = [0, py, 0]
    jacobian[4] = {0.0, py, 0.0};
    // dy/dM12 = [0, pz, 0]
    jacobian[5] = {0.0, pz, 0.0};
    
    // dy/dM20 = [0, 0, px]
    jacobian[6] = {0.0, 0.0, px};
    // dy/dM21 = [0, 0, py]
    jacobian[7] = {0.0, 0.0, py};
    // dy/dM22 = [0, 0, pz]
    jacobian[8] = {0.0, 0.0, pz};
    
    // dy/dTx = [1, 0, 0]
    jacobian[9] = {1.0, 0.0, 0.0};
    // dy/dTy = [0, 1, 0]
    jacobian[10] = {0.0, 1.0, 0.0};
    // dy/dTz = [0, 0, 1]
    jacobian[11] = {0.0, 0.0, 1.0};
}

// ============================================================================
// 单层配准 - 刚体
// ============================================================================

void ImageRegistration::RunSingleLevelRigid(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level)
{
    // 获取当前层的迭代次数
    unsigned int currentIterations = m_NumberOfIterations[0];  // 默认使用第一个值
    if (level < m_NumberOfIterations.size())
    {
        currentIterations = m_NumberOfIterations[level];
    }
    
    // 如果当前层迭代次数为0,直接跳过
    if (currentIterations == 0)
    {
        std::cout << "  [Skipping] Level " << level << " iterations set to 0" << std::endl;
        return;
    }
    
    // 配置度量
    m_Metric->SetFixedImage(fixedImage);
    m_Metric->SetMovingImage(movingImage);
    m_Metric->SetNumberOfHistogramBins(m_NumberOfHistogramBins);
    if (m_SamplingPercentage > 0.0 && m_SamplingPercentage <= 1.0)
    {
        m_Metric->SetSamplingPercentage(m_SamplingPercentage);
    }
    else
    {
        m_Metric->SetNumberOfSpatialSamples(m_NumberOfSpatialSamples);
    }
    m_Metric->SetRandomSeed(m_RandomSeed);
    
    // 设置掩膜 (如果有的话,用于局部配准)
    if (m_FixedImageMask.IsNotNull())
    {
        m_Metric->SetFixedImageMask(m_FixedImageMask);
    }
    
    // 关键修复：直接使用已经从初始变换初始化好参数的 m_RigidTransform
    // 不再使用 CompositeTransform！这是 ANTs 的做法
    m_Metric->SetTransform(m_RigidTransform);
    m_Metric->SetNumberOfParameters(6);
    m_Metric->SetUseStratifiedSampling(m_UseStratifiedSampling);
    
    // 设置雅可比函数(不需要链式法则，因为没有 composite)
    auto rigidTransformPtr = m_RigidTransform;
    m_Metric->SetJacobianFunction([rigidTransformPtr](const ImageType::PointType& point,
                                                       std::vector<std::array<double, 3>>& jacobian) {
        ComputeRigidJacobian(point, rigidTransformPtr, jacobian);
    });
    
    m_Metric->Initialize();

    // 配置优化器 - 使用分层学习率
    double currentLearningRate = (level < m_LearningRate.size()) 
                                  ? m_LearningRate[level] 
                                  : m_LearningRate.back();
    std::cout << "  Learning Rate: " << std::fixed << std::setprecision(4) << currentLearningRate << std::endl;
    
    m_Optimizer->SetLearningRate(currentLearningRate);
    m_Optimizer->SetMinimumStepLength(m_MinimumStepLength);
    m_Optimizer->SetNumberOfIterations(currentIterations);  // 使用当前层的迭代次数
    m_Optimizer->SetRelaxationFactor(m_RelaxationFactor);
    m_Optimizer->SetGradientMagnitudeTolerance(m_GradientMagnitudeTolerance);
    m_Optimizer->SetReturnBestParametersAndValue(true);
    m_Optimizer->SetNumberOfParameters(6);
    
    std::vector<double> scales = EstimateRigidParameterScales();
    m_Optimizer->SetScales(scales);

    m_Optimizer->SetCostFunction([this]() -> double {
        return m_Metric->GetValue();
    });

    m_Optimizer->SetGradientFunction([this](std::vector<double>& gradient) {
        m_Metric->GetDerivative(gradient);
    });
    
    m_Optimizer->SetGetParametersFunction([this]() -> std::vector<double> {
        auto params = m_RigidTransform->GetParameters();
        std::vector<double> result(params.Size());
        for (unsigned int i = 0; i < params.Size(); ++i)
        {
            result[i] = params[i];
        }
        return result;
    });
    
    m_Optimizer->SetSetParametersFunction([this](const std::vector<double>& params) {
        RigidTransformType::ParametersType itkParams(6);
        for (unsigned int i = 0; i < 6 && i < params.size(); ++i)
        {
            itkParams[i] = params[i];
        }
        m_RigidTransform->SetParameters(itkParams);
    });

    m_Optimizer->SetUpdateParametersFunction([this](const std::vector<double>& update) {
        auto params = m_RigidTransform->GetParameters();
        for (unsigned int i = 0; i < 6 && i < params.Size(); ++i)
        {
            params[i] += update[i];
        }
        m_RigidTransform->SetParameters(params);
    });

    // 每次迭代打印观察者间隔(如 verbose 则打印每次迭代)
    m_Optimizer->SetVerbose(m_Verbose);
    if (m_Verbose)
    {
        m_Optimizer->SetObserverIterationInterval(1);
    }
    else
    {
        m_Optimizer->SetObserverIterationInterval(10);
    }

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

    m_Optimizer->StartOptimization();
    m_FinalMetricValue = m_Optimizer->GetBestValue();
}

// ============================================================================
// 单层配准 - 仿射
// ============================================================================

void ImageRegistration::RunSingleLevelAffine(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level)
{
    // 获取当前层的迭代次数
    unsigned int currentIterations = m_NumberOfIterations[0];
    if (level < m_NumberOfIterations.size())
    {
        currentIterations = m_NumberOfIterations[level];
    }
    
    // 如果当前层迭代次数为0,直接跳过
    if (currentIterations == 0)
    {
        std::cout << "  [Skipping] Level " << level << " iterations set to 0" << std::endl;
        return;
    }
    
    // 配置度量
    m_Metric->SetFixedImage(fixedImage);
    m_Metric->SetMovingImage(movingImage);
    m_Metric->SetNumberOfHistogramBins(m_NumberOfHistogramBins);
    if (m_SamplingPercentage > 0.0 && m_SamplingPercentage <= 1.0)
    {
        m_Metric->SetSamplingPercentage(m_SamplingPercentage);
    }
    else
    {
        m_Metric->SetNumberOfSpatialSamples(m_NumberOfSpatialSamples);
    }
    m_Metric->SetRandomSeed(m_RandomSeed);
    
    // 设置掩膜 (如果有的话,用于局部配准)
    if (m_FixedImageMask.IsNotNull())
    {
        m_Metric->SetFixedImageMask(m_FixedImageMask);
    }
    
    // 关键:直接使用已初始化的m_AffineTransform,不使用CompositeTransform
    // 与刚体配准采用相同的策略
    m_Metric->SetTransform(m_AffineTransform);
    m_Metric->SetNumberOfParameters(12);
    m_Metric->SetUseStratifiedSampling(m_UseStratifiedSampling);
    
    // 设置雅可比函数
    auto affineTransformPtr = m_AffineTransform;
    m_Metric->SetJacobianFunction([affineTransformPtr](const ImageType::PointType& point,
                                                       std::vector<std::array<double, 3>>& jacobian) {
        ComputeAffineJacobian(point, affineTransformPtr, jacobian);
    });
    
    m_Metric->Initialize();

    // 配置优化器 - 使用分层学习率
    double currentLearningRate = (level < m_LearningRate.size()) 
                                  ? m_LearningRate[level] 
                                  : m_LearningRate.back();
    std::cout << "  Learning Rate: " << std::fixed << std::setprecision(4) << currentLearningRate << std::endl;
    
    m_Optimizer->SetLearningRate(currentLearningRate);
    m_Optimizer->SetMinimumStepLength(m_MinimumStepLength);
    m_Optimizer->SetNumberOfIterations(currentIterations);  // 使用当前层的迭代次数
    m_Optimizer->SetRelaxationFactor(m_RelaxationFactor);
    m_Optimizer->SetGradientMagnitudeTolerance(m_GradientMagnitudeTolerance);
    m_Optimizer->SetReturnBestParametersAndValue(true);
    m_Optimizer->SetNumberOfParameters(12);
    
    std::vector<double> scales = EstimateAffineParameterScales();
    m_Optimizer->SetScales(scales);
    
    // 为仿射变换设置合理的参数更新上限
    // 矩阵元素:每次迭代最多变化0.1 (10%形变)
    // 平移:每次迭代最多20mm
    std::vector<double> maxUpdate(12);
    for (int i = 0; i < 9; ++i)
    {
        maxUpdate[i] = 0.1;  // 矩阵元素
    }
    for (int i = 9; i < 12; ++i)
    {
        maxUpdate[i] = 20.0;  // 平移(mm)
    }
    m_Optimizer->SetMaxParameterUpdate(maxUpdate);

    m_Optimizer->SetCostFunction([this]() -> double {
        return m_Metric->GetValue();
    });

    m_Optimizer->SetGradientFunction([this](std::vector<double>& gradient) {
        m_Metric->GetDerivative(gradient);
    });
    
    m_Optimizer->SetGetParametersFunction([this]() -> std::vector<double> {
        auto params = m_AffineTransform->GetParameters();
        std::vector<double> result(params.Size());
        for (unsigned int i = 0; i < params.Size(); ++i)
        {
            result[i] = params[i];
        }
        return result;
    });
    
    m_Optimizer->SetSetParametersFunction([this](const std::vector<double>& params) {
        AffineTransformType::ParametersType itkParams(12);
        for (unsigned int i = 0; i < 12 && i < params.size(); ++i)
        {
            itkParams[i] = params[i];
        }
        m_AffineTransform->SetParameters(itkParams);
    });

    m_Optimizer->SetUpdateParametersFunction([this](const std::vector<double>& update) {
        auto params = m_AffineTransform->GetParameters();
        for (unsigned int i = 0; i < 12 && i < params.Size(); ++i)
        {
            params[i] += update[i];
        }
        m_AffineTransform->SetParameters(params);
    });

    // 设置观察者
    m_Optimizer->SetVerbose(m_Verbose);
    if (m_Verbose)
    {
        m_Optimizer->SetObserverIterationInterval(1);
    }
    else
    {
        m_Optimizer->SetObserverIterationInterval(10);
    }

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

    m_Optimizer->StartOptimization();
    m_FinalMetricValue = m_Optimizer->GetBestValue();
}

// ============================================================================
// 主配准流程
// ============================================================================

void ImageRegistration::Update()
{
    if (!m_FixedImage || !m_MovingImage)
    {
        throw std::runtime_error("Fixed or moving image not set");
    }

    auto startTime = std::chrono::high_resolution_clock::now();

    // 初始化变换
    InitializeTransform();

    // 把 verbose 传递下去
    m_Metric->SetVerbose(m_Verbose);
    m_Optimizer->SetVerbose(m_Verbose);

    // 打印变换类型
    std::cout << "\nTransform Type: " << ConfigManager::TransformTypeToString(m_TransformType) 
              << " (" << GetNumberOfParameters() << " parameters)" << std::endl;

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

    // 计算采样信息 (考虑掩膜)
    auto fixedRegion = m_FixedImage->GetLargestPossibleRegion();
    auto fixedSize = fixedRegion.GetSize();
    unsigned long totalVoxels = fixedSize[0] * fixedSize[1] * fixedSize[2];
    unsigned long effectiveVoxels = totalVoxels;  // 有效体素数(考虑掩膜)
    
    // 如果有掩膜,使用之前统计好的掩膜体素数
    if (m_FixedImageMask.IsNotNull() && m_MaskVoxelCount > 0)
    {
        effectiveVoxels = m_MaskVoxelCount;
    }
    
    unsigned int effectiveSamples = m_NumberOfSpatialSamples;
    double samplingPctDisplay = m_SamplingPercentage * 100.0;
    if (m_SamplingPercentage > 0.0 && m_SamplingPercentage <= 1.0)
    {
        unsigned long computed = static_cast<unsigned long>(std::round(m_SamplingPercentage * static_cast<double>(effectiveVoxels)));
        if (computed < 1) computed = 1;
        if (computed > effectiveVoxels) computed = effectiveVoxels;
        effectiveSamples = static_cast<unsigned int>(computed);
        samplingPctDisplay = m_SamplingPercentage * 100.0;
    }
    else
    {
        // compute percent from explicit samples
        samplingPctDisplay = static_cast<double>(effectiveSamples) / static_cast<double>(effectiveVoxels) * 100.0;
    }
    
    if (m_FixedImageMask.IsNotNull())
    {
        std::cout << "Using " << effectiveSamples << " samples (" 
                  << std::fixed << std::setprecision(1) << samplingPctDisplay 
                  << "% of mask region, " << effectiveVoxels << " voxels in mask)\n";
    }
    else
    {
        std::cout << "Using " << effectiveSamples << " samples (" 
                  << std::fixed << std::setprecision(1) << samplingPctDisplay 
                  << "% of total voxels)\n";
    }

    // 多分辨率金字塔
    for (unsigned int level = 0; level < m_NumberOfLevels; ++level)
    {
        unsigned int shrinkFactor = (level < m_ShrinkFactors.size()) ? m_ShrinkFactors[level] : 1;
        double smoothingSigma = (level < m_SmoothingSigmas.size()) ? m_SmoothingSigmas[level] : 0.0;
        
        if (m_LevelObserver)
        {
            m_LevelObserver(level, shrinkFactor, smoothingSigma);
        }
        else
        {
            std::cout << "\n[Multi-Resolution Level " << (level + 1) << "]" << std::endl;
            std::cout << "  Shrink Factor: " << shrinkFactor << "x" << std::endl;
            std::cout << "  Smoothing Sigma: " << std::fixed << std::setprecision(2) << smoothingSigma << " mm" << std::endl;
        }

        // ANTs风格的预处理：Winsorizing -> Smooth -> Shrink
        ImageType::Pointer fixedPyramid = WinsorizeImage(m_FixedImage, 0.005, 0.995);
        fixedPyramid = SmoothImage(fixedPyramid, smoothingSigma);
        fixedPyramid = ShrinkImage(fixedPyramid, shrinkFactor);

        ImageType::Pointer movingPyramid = WinsorizeImage(m_MovingImage, 0.005, 0.995);
        movingPyramid = SmoothImage(movingPyramid, smoothingSigma);
        movingPyramid = ShrinkImage(movingPyramid, shrinkFactor);

        RunSingleLevel(fixedPyramid, movingPyramid, level);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    m_ElapsedTime = std::chrono::duration<double>(endTime - startTime).count();

    std::cout << "Final metric value: " << std::scientific << std::setprecision(4) 
              << m_FinalMetricValue << std::endl;
}

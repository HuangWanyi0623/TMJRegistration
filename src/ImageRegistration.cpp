#include "ImageRegistration.h"
#include "itkIdentityTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <iostream>

ImageRegistration::ImageRegistration()
    : m_NumberOfHistogramBins(50)
    , m_NumberOfSpatialSamples(10000)
    , m_MaximumStepLength(1.0)
    , m_MinimumStepLength(0.001)
    , m_NumberOfIterations(200)
    , m_FinalMetricValue(0.0)
{
    m_Metric = std::make_unique<MattesMutualInformation>();
    m_Optimizer = std::make_unique<RegularStepGradientDescentOptimizer>();
}

ImageRegistration::~ImageRegistration()
{
}

void ImageRegistration::SetFixedImage(ImageType::Pointer fixedImage)
{
    m_FixedImage = fixedImage;
}

void ImageRegistration::SetMovingImage(ImageType::Pointer movingImage)
{
    m_MovingImage = movingImage;
}

void ImageRegistration::InitializeTransformParameters()
{
    // 初始化为6参数刚体变换 (3平移 + 3旋转)
    m_FinalParameters.resize(6, 0.0);
    
    // 可以在这里设置初始估计值
    // 目前设置为单位变换(全0)
}

void ImageRegistration::Update()
{
    if (!m_FixedImage || !m_MovingImage)
    {
        throw std::runtime_error("Fixed or moving image not set");
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "   Image Registration Started" << std::endl;
    std::cout << "========================================" << std::endl;

    // 打印图像信息
    auto fixedRegion = m_FixedImage->GetLargestPossibleRegion();
    auto movingRegion = m_MovingImage->GetLargestPossibleRegion();
    
    std::cout << "\nFixed Image:" << std::endl;
    std::cout << "  Size: " << fixedRegion.GetSize() << std::endl;
    std::cout << "  Spacing: " << m_FixedImage->GetSpacing() << std::endl;
    std::cout << "  Origin: " << m_FixedImage->GetOrigin() << std::endl;

    std::cout << "\nMoving Image:" << std::endl;
    std::cout << "  Size: " << movingRegion.GetSize() << std::endl;
    std::cout << "  Spacing: " << m_MovingImage->GetSpacing() << std::endl;
    std::cout << "  Origin: " << m_MovingImage->GetOrigin() << std::endl;

    // 配置度量
    m_Metric->SetFixedImage(m_FixedImage);
    m_Metric->SetMovingImage(m_MovingImage);
    m_Metric->SetNumberOfHistogramBins(m_NumberOfHistogramBins);
    m_Metric->SetNumberOfSpatialSamples(m_NumberOfSpatialSamples);
    
    std::cout << "\n--- Initializing Metric ---" << std::endl;
    m_Metric->Initialize();

    // 初始化变换参数
    InitializeTransformParameters();

    // 配置优化器
    m_Optimizer->SetMaximumStepLength(m_MaximumStepLength);
    m_Optimizer->SetMinimumStepLength(m_MinimumStepLength);
    m_Optimizer->SetNumberOfIterations(m_NumberOfIterations);
    m_Optimizer->SetRelaxationFactor(0.5);
    m_Optimizer->SetGradientMagnitudeTolerance(1e-4);
    m_Optimizer->SetInitialPosition(m_FinalParameters);

    // 设置代价函数
    auto costFunction = [this](const ParametersType& params) -> double {
        return m_Metric->GetValue(params);
    };
    m_Optimizer->SetCostFunction(costFunction);

    // 设置梯度函数
    auto gradientFunction = [this](const ParametersType& params, ParametersType& gradient) {
        m_Metric->GetDerivative(params, gradient);
    };
    m_Optimizer->SetGradientFunction(gradientFunction);

    // 执行优化
    std::cout << "\n--- Starting Optimization ---" << std::endl;
    m_Optimizer->StartOptimization();

    // 获取最终结果
    m_FinalParameters = m_Optimizer->GetCurrentPosition();
    m_FinalMetricValue = m_Optimizer->GetCurrentValue();

    std::cout << "\n========================================" << std::endl;
    std::cout << "   Registration Results" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Final Parameters:" << std::endl;
    std::cout << "  Translation: [" << m_FinalParameters[0] << ", " 
              << m_FinalParameters[1] << ", " << m_FinalParameters[2] << "]" << std::endl;
    std::cout << "  Rotation (rad): [" << m_FinalParameters[3] << ", " 
              << m_FinalParameters[4] << ", " << m_FinalParameters[5] << "]" << std::endl;
    std::cout << "  Rotation (deg): [" << m_FinalParameters[3] * 180.0 / 3.14159 << ", " 
              << m_FinalParameters[4] * 180.0 / 3.14159 << ", " 
              << m_FinalParameters[5] * 180.0 / 3.14159 << "]" << std::endl;
    std::cout << "Final Metric Value: " << m_FinalMetricValue << std::endl;

    // 应用变换生成输出图像
    std::cout << "\n--- Generating Output Image ---" << std::endl;
    
    using ResampleFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    auto resampler = ResampleFilterType::New();

    // 创建仿射变换
    auto transform = TransformType::New();
    TransformType::ParametersType transformParams(12);
    
    // 提取旋转和平移
    double tx = m_FinalParameters[0];
    double ty = m_FinalParameters[1];
    double tz = m_FinalParameters[2];
    double rx = m_FinalParameters[3];
    double ry = m_FinalParameters[4];
    double rz = m_FinalParameters[5];

    // 计算旋转矩阵
    double cx = std::cos(rx), sx = std::sin(rx);
    double cy = std::cos(ry), sy = std::sin(ry);
    double cz = std::cos(rz), sz = std::sin(rz);

    // 设置变换矩阵参数(ITK使用逆变换)
    transformParams[0] = cy * cz;
    transformParams[1] = -cy * sz;
    transformParams[2] = sy;
    transformParams[3] = sx * sy * cz + cx * sz;
    transformParams[4] = -sx * sy * sz + cx * cz;
    transformParams[5] = -sx * cy;
    transformParams[6] = -cx * sy * cz + sx * sz;
    transformParams[7] = cx * sy * sz + sx * cz;
    transformParams[8] = cx * cy;
    transformParams[9] = tx;
    transformParams[10] = ty;
    transformParams[11] = tz;

    transform->SetParameters(transformParams);

    // 配置重采样滤波器
    resampler->SetTransform(transform);
    resampler->SetInput(m_MovingImage);
    resampler->SetSize(m_FixedImage->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(m_FixedImage->GetOrigin());
    resampler->SetOutputSpacing(m_FixedImage->GetSpacing());
    resampler->SetOutputDirection(m_FixedImage->GetDirection());
    resampler->SetDefaultPixelValue(0);

    try
    {
        resampler->Update();
        m_OutputImage = resampler->GetOutput();
        std::cout << "Output image generated successfully" << std::endl;
    }
    catch (const itk::ExceptionObject& e)
    {
        std::cerr << "Error during resampling: " << e << std::endl;
        throw;
    }
}

ImageRegistration::ImageType::Pointer ImageRegistration::GetOutput()
{
    return m_OutputImage;
}

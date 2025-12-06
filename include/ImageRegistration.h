#ifndef IMAGE_REGISTRATION_H
#define IMAGE_REGISTRATION_H

#include "MattesMutualInformation.h"
#include "RegularStepGradientDescentOptimizer.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"

/**
 * @brief 图像配准主类
 * 
 * 整合互信息度量和梯度下降优化器,实现3D图像配准
 */
class ImageRegistration
{
public:
    using ImageType = itk::Image<float, 3>;
    using TransformType = itk::AffineTransform<double, 3>;
    using ParametersType = std::vector<double>;

    ImageRegistration();
    ~ImageRegistration();

    // 设置图像
    void SetFixedImage(ImageType::Pointer fixedImage);
    void SetMovingImage(ImageType::Pointer movingImage);

    // 设置配准参数
    void SetNumberOfHistogramBins(unsigned int bins) { m_NumberOfHistogramBins = bins; }
    void SetNumberOfSpatialSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; }
    void SetMaximumStepLength(double stepLength) { m_MaximumStepLength = stepLength; }
    void SetMinimumStepLength(double stepLength) { m_MinimumStepLength = stepLength; }
    void SetNumberOfIterations(unsigned int iterations) { m_NumberOfIterations = iterations; }

    // 执行配准
    void Update();

    // 获取配准后的图像
    ImageType::Pointer GetOutput();

    // 获取变换参数
    const ParametersType& GetTransformParameters() const { return m_FinalParameters; }

    // 获取最终度量值
    double GetFinalMetricValue() const { return m_FinalMetricValue; }

private:
    // 输入图像
    ImageType::Pointer m_FixedImage;
    ImageType::Pointer m_MovingImage;
    ImageType::Pointer m_OutputImage;

    // 度量和优化器
    std::unique_ptr<MattesMutualInformation> m_Metric;
    std::unique_ptr<RegularStepGradientDescentOptimizer> m_Optimizer;

    // 配准参数
    unsigned int m_NumberOfHistogramBins;
    unsigned int m_NumberOfSpatialSamples;
    double m_MaximumStepLength;
    double m_MinimumStepLength;
    unsigned int m_NumberOfIterations;

    // 结果
    ParametersType m_FinalParameters;
    double m_FinalMetricValue;

    // 初始化变换参数
    void InitializeTransformParameters();
};

#endif // IMAGE_REGISTRATION_H

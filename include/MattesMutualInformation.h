#ifndef MATTES_MUTUAL_INFORMATION_H
#define MATTES_MUTUAL_INFORMATION_H

#include <vector>
#include <memory>
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"

/**
 * @brief Mattes互信息度量类
 * 
 * 基于ITK的itkMattesMutualInformationImageToImageMetric实现
 * 使用Parzen窗方法估计联合概率分布和边缘概率分布
 */
class MattesMutualInformation
{
public:
    using ImageType = itk::Image<float, 3>;
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    using ParametersType = std::vector<double>;

    MattesMutualInformation();
    ~MattesMutualInformation();

    // 设置固定图像和移动图像
    void SetFixedImage(ImageType::Pointer fixedImage);
    void SetMovingImage(ImageType::Pointer movingImage);

    // 设置参数
    void SetNumberOfHistogramBins(unsigned int bins) { m_NumberOfHistogramBins = bins; }
    void SetNumberOfSpatialSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; }

    // 初始化
    void Initialize();

    // 计算互信息值和梯度
    double GetValue(const ParametersType& parameters);
    void GetDerivative(const ParametersType& parameters, ParametersType& derivative);
    void GetValueAndDerivative(const ParametersType& parameters, double& value, ParametersType& derivative);

private:
    // 图像指针
    ImageType::Pointer m_FixedImage;
    ImageType::Pointer m_MovingImage;
    InterpolatorType::Pointer m_Interpolator;

    // 直方图参数
    unsigned int m_NumberOfHistogramBins;
    unsigned int m_NumberOfSpatialSamples;

    // 联合直方图和边缘直方图
    std::vector<std::vector<double>> m_JointPDF;
    std::vector<double> m_FixedImageMarginalPDF;
    std::vector<double> m_MovingImageMarginalPDF;

    // 采样点
    struct SamplePoint
    {
        ImageType::PointType point;
        double fixedValue;
    };
    std::vector<SamplePoint> m_SamplePoints;

    // 图像强度范围
    double m_FixedImageMin;
    double m_FixedImageMax;
    double m_MovingImageMin;
    double m_MovingImageMax;

    // 内部方法
    void ComputeImageExtrema();
    void SampleFixedImage();
    void ComputeJointPDF(const ParametersType& parameters);
    double ComputeMutualInformation();
    void ApplyTransform(const ImageType::PointType& inputPoint, 
                       const ParametersType& parameters,
                       ImageType::PointType& outputPoint);
    int ComputeFixedImageBin(double value) const;
    int ComputeMovingImageBin(double value) const;
};

#endif // MATTES_MUTUAL_INFORMATION_H

#ifndef MATTES_MUTUAL_INFORMATION_H
#define MATTES_MUTUAL_INFORMATION_H

#include <vector>
#include <memory>
#include <random>
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"

/**
 * @brief Mattes互信息度量类
 * 
 * 基于ITK的itkMattesMutualInformationImageToImageMetric实现
 * 使用Parzen窗方法估计联合概率分布和边缘概率分布
 * 
 * 扩展性说明:
 * - 可以修改ComputeJointPDF()改变直方图计算方式
 * - 可以修改ComputeMutualInformation()使用不同的信息度量
 * - 可以添加B样条Parzen窗估计替代简单直方图
 */
class MattesMutualInformation
{
public:
    using ImageType = itk::Image<float, 3>;
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    using TransformType = itk::Euler3DTransform<double>;
    using ParametersType = std::vector<double>;

    MattesMutualInformation();
    ~MattesMutualInformation();

    // 设置固定图像和移动图像
    void SetFixedImage(ImageType::Pointer fixedImage);
    void SetMovingImage(ImageType::Pointer movingImage);
    
    // 设置变换(用于物理空间计算)
    void SetTransform(TransformType::Pointer transform) { m_Transform = transform; }

    // 设置参数
    void SetNumberOfHistogramBins(unsigned int bins) { m_NumberOfHistogramBins = bins; }
    void SetNumberOfSpatialSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; }
    void SetRandomSeed(unsigned int seed) { m_RandomSeed = seed; m_UseFixedSeed = true; }

    // 初始化
    void Initialize();
    
    // 重新采样(用于多分辨率)
    void ReinitializeSampling();

    // 计算互信息值和梯度
    double GetValue();
    void GetDerivative(ParametersType& derivative);
    void GetValueAndDerivative(double& value, ParametersType& derivative);

    // 获取当前度量值
    double GetCurrentValue() const { return m_CurrentValue; }

private:
    // 图像指针
    ImageType::Pointer m_FixedImage;
    ImageType::Pointer m_MovingImage;
    InterpolatorType::Pointer m_Interpolator;
    TransformType::Pointer m_Transform;

    // 直方图参数
    unsigned int m_NumberOfHistogramBins;
    unsigned int m_NumberOfSpatialSamples;
    unsigned int m_RandomSeed;
    bool m_UseFixedSeed;

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

    // 当前度量值
    double m_CurrentValue;

    // 随机数生成器
    std::mt19937 m_RandomGenerator;

    // 内部方法
    void ComputeImageExtrema();
    void SampleFixedImage();
    void ComputeJointPDF();
    double ComputeMutualInformation();
    int ComputeFixedImageBin(double value) const;
    int ComputeMovingImageBin(double value) const;
};

#endif // MATTES_MUTUAL_INFORMATION_H

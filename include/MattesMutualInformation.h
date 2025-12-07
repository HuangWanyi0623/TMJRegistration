#ifndef MATTES_MUTUAL_INFORMATION_H
#define MATTES_MUTUAL_INFORMATION_H

#include <vector>
#include <memory>
#include <random>
#include <array>
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"

/**
 * @brief Mattes互信息度量类 - 带解析梯度
 * 
 * 基于ITK的itkMattesMutualInformationImageToImageMetricv4实现
 * 核心特性:
 * - 使用三次B样条Parzen窗估计概率密度
 * - 实现解析梯度计算(非有限差分)
 * - 支持均匀分层采样策略
 * 
 * 数学原理:
 * MI = H(F) + H(M) - H(F,M)
 * dMI/dp = -sum_f sum_m [ (1 + log(p(f,m)/p(m))) * dp(f,m)/dp ]
 * 
 * 其中 dp(f,m)/dp 通过链式法则计算:
 * dp(f,m)/dp = sum_samples [ dB/dm * dm/dp ]
 * dm/dp = gradient_M(T(x)) * dT/dp
 */
class MattesMutualInformation
{
public:
    using ImageType = itk::Image<float, 3>;
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    using TransformType = itk::Euler3DTransform<double>;
    using ParametersType = std::vector<double>;
    
    // B样条阶数 (ITK使用3阶)
    static constexpr unsigned int BSplineOrder = 3;
    static constexpr unsigned int NumberOfBSplineCoefficients = BSplineOrder + 1; // 4

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
    
    // 采样策略设置
    void SetUseStratifiedSampling(bool use) { m_UseStratifiedSampling = use; }

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
    
    // 获取有效采样点数量(用于调试)
    unsigned int GetNumberOfValidSamples() const { return m_NumberOfValidSamples; }

private:
    // 图像指针
    ImageType::Pointer m_FixedImage;
    ImageType::Pointer m_MovingImage;
    InterpolatorType::Pointer m_Interpolator;
    TransformType::Pointer m_Transform;

    // 移动图像梯度(用于解析梯度计算)
    std::array<ImageType::Pointer, 3> m_MovingImageGradient;
    std::array<InterpolatorType::Pointer, 3> m_GradientInterpolators;

    // 直方图参数
    unsigned int m_NumberOfHistogramBins;
    unsigned int m_NumberOfSpatialSamples;
    unsigned int m_RandomSeed;
    bool m_UseFixedSeed;
    bool m_UseStratifiedSampling;
    unsigned int m_NumberOfValidSamples;

    // B样条Parzen窗的联合直方图 (二维)
    std::vector<std::vector<double>> m_JointPDF;
    
    // 边缘概率分布
    std::vector<double> m_FixedImageMarginalPDF;
    std::vector<double> m_MovingImageMarginalPDF;
    
    // 用于解析梯度的联合PDF导数存储
    // m_JointPDFDerivatives[param][fixedBin][movingBin]
    std::vector<std::vector<std::vector<double>>> m_JointPDFDerivatives;

    // 采样点信息(扩展版)
    struct SamplePoint
    {
        ImageType::PointType fixedPoint;     // 固定图像中的物理点
        double fixedValue;                   // 固定图像强度
        int fixedParzenWindowIndex;          // 固定图像B样条起始索引
        std::array<double, 4> fixedBSplineWeights;  // 固定图像B样条权重
    };
    std::vector<SamplePoint> m_SamplePoints;

    // 图像强度范围
    double m_FixedImageMin;
    double m_FixedImageMax;
    double m_MovingImageMin;
    double m_MovingImageMax;
    
    // 强度到bin的转换系数
    double m_FixedImageBinSize;
    double m_MovingImageBinSize;

    // 当前度量值
    double m_CurrentValue;

    // 随机数生成器
    std::mt19937 m_RandomGenerator;

    // ============ 内部方法 ============
    
    // 初始化相关
    void ComputeImageExtrema();
    void ComputeMovingImageGradient();
    
    // 采样策略
    void SampleFixedImage();
    void SampleFixedImageStratified();  // 分层均匀采样
    void SampleFixedImageRandom();      // 随机采样
    
    // B样条相关
    double EvaluateCubicBSpline(double u) const;
    double EvaluateCubicBSplineDerivative(double u) const;
    void ComputeBSplineWeights(double continuousIndex, int& startIndex, 
                               std::array<double, 4>& weights) const;
    void ComputeBSplineDerivativeWeights(double continuousIndex, int& startIndex,
                                         std::array<double, 4>& derivativeWeights) const;
    
    // 核心计算
    void ComputeJointPDFAndDerivatives();
    double ComputeMutualInformation();
    void ComputeAnalyticalGradient(ParametersType& derivative);
    
    // 变换雅可比矩阵
    void ComputeTransformJacobian(const ImageType::PointType& point,
                                  std::vector<std::array<double, 3>>& jacobian);
    
    // 辅助函数
    double ComputeFixedImageContinuousIndex(double value) const;
    double ComputeMovingImageContinuousIndex(double value) const;
};

#endif // MATTES_MUTUAL_INFORMATION_H

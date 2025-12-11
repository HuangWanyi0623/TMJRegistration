#ifndef MATTES_MUTUAL_INFORMATION_H
#define MATTES_MUTUAL_INFORMATION_H

#include <vector>
#include <memory>
#include <random>
#include <array>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTransform.h"
#include "itkImageMaskSpatialObject.h"

/**
 * @brief Mattes互信息度量类 - 带解析梯度
 * 
 * 基于ITK的itkMattesMutualInformationImageToImageMetricv4实现
 * 核心特性:
 * - 使用三次B样条Parzen窗估计概率密度
 * - 实现解析梯度计算(非有限差分)
 * - 支持均匀分层采样策略
 * - 支持任意维度的变换参数(刚体6参数/仿射12参数)
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
    using MaskImageType = itk::Image<unsigned char, 3>;
    using MaskSpatialObjectType = itk::ImageMaskSpatialObject<3>;
    using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
    using TransformBaseType = itk::Transform<double, 3, 3>;
    using ParametersType = std::vector<double>;
    
    // 雅可比矩阵计算回调类型
    using JacobianFunctionType = std::function<void(const ImageType::PointType&, 
                                                     std::vector<std::array<double, 3>>&)>;
    
    // B样条阶数 (ITK使用3阶)
    static constexpr unsigned int BSplineOrder = 3;
    static constexpr unsigned int NumberOfBSplineCoefficients = BSplineOrder + 1; // 4

    MattesMutualInformation();
    ~MattesMutualInformation();

    void SetVerbose(bool v) { m_Verbose = v; }
    bool GetVerbose() const { return m_Verbose; }

    // 设置固定图像和移动图像
    void SetFixedImage(ImageType::Pointer fixedImage);
    void SetMovingImage(ImageType::Pointer movingImage);
    
    // 设置变换(使用通用变换基类)
    void SetTransform(TransformBaseType::Pointer transform) { m_Transform = transform; }
    
    // 设置雅可比矩阵计算函数(由外部提供,支持不同变换类型)
    void SetJacobianFunction(JacobianFunctionType func) { m_JacobianFunction = func; }
    
    // 设置参数数量(根据变换类型: 刚体6, 仿射12)
    void SetNumberOfParameters(unsigned int num) { m_NumberOfParameters = num; }

    // 设置参数
    void SetNumberOfHistogramBins(unsigned int bins) { m_NumberOfHistogramBins = bins; }
    void SetNumberOfSpatialSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; m_SamplingPercentage = 0.0; }
    // 采样百分比 [0.0, 1.0], 当>0时优先使用百分比计算采样数 (默认 0.10)
    void SetSamplingPercentage(double percent) { m_SamplingPercentage = percent; if (percent > 0.0) m_NumberOfSpatialSamples = 0; }
    double GetSamplingPercentage() const { return m_SamplingPercentage; }
    void SetRandomSeed(unsigned int seed) { m_RandomSeed = seed; m_UseFixedSeed = true; }
    
    // 掩膜设置 (用于局部配准,只在掩膜区域内采样)
    void SetFixedImageMask(MaskSpatialObjectType::Pointer mask) { m_FixedImageMask = mask; }
    MaskSpatialObjectType::Pointer GetFixedImageMask() const { return m_FixedImageMask; }
    bool HasFixedImageMask() const { return m_FixedImageMask.IsNotNull(); }
    
    // 采样策略设置
    void SetUseStratifiedSampling(bool use) { m_UseStratifiedSampling = use; }
    
    // 多线程设置
    void SetNumberOfThreads(unsigned int n) { m_NumberOfThreads = n; }
    unsigned int GetNumberOfThreads() const { return m_NumberOfThreads; }

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
    TransformBaseType::Pointer m_Transform;
    
    // 掩膜 (可选,用于局部配准)
    MaskSpatialObjectType::Pointer m_FixedImageMask;
    
    // 雅可比矩阵计算函数(外部提供)
    JacobianFunctionType m_JacobianFunction;
    unsigned int m_NumberOfParameters;

    // 移动图像梯度(用于解析梯度计算)
    std::array<ImageType::Pointer, 3> m_MovingImageGradient;
    std::array<InterpolatorType::Pointer, 3> m_GradientInterpolators;

    // 直方图参数
    unsigned int m_NumberOfHistogramBins;
    unsigned int m_NumberOfSpatialSamples;
    double m_SamplingPercentage;
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
    bool m_Verbose; // 调试输出
    
    // 多线程参数
    unsigned int m_NumberOfThreads;
    
    // 多线程局部直方图 (每个线程一个)
    struct ThreadLocalHistograms
    {
        std::vector<std::vector<double>> jointPDF;
        std::vector<std::vector<std::vector<double>>> jointPDFDerivatives;
        unsigned int validSamples;
        
        ThreadLocalHistograms(unsigned int numBins, unsigned int numParams)
            : validSamples(0)
        {
            jointPDF.resize(numBins, std::vector<double>(numBins, 0.0));
            jointPDFDerivatives.resize(numParams);
            for (auto& paramDeriv : jointPDFDerivatives)
            {
                paramDeriv.resize(numBins, std::vector<double>(numBins, 0.0));
            }
        }
    };

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
    void ComputeJointPDFAndDerivativesThreaded();  // 多线程版本
    void ComputePDFRange(size_t startIdx, size_t endIdx, ThreadLocalHistograms& localHist);
    double ComputeMutualInformation();
    void ComputeAnalyticalGradient(ParametersType& derivative);
    
    // 辅助函数
    double ComputeFixedImageContinuousIndex(double value) const;
    double ComputeMovingImageContinuousIndex(double value) const;
};

#endif // MATTES_MUTUAL_INFORMATION_H

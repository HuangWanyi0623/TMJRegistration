#ifndef IMAGEREGISTRATION_H
#define IMAGEREGISTRATION_H

#include <vector>
#include <memory>
#include <functional>
#include <string>
#include <itkImage.h>
#include <itkEuler3DTransform.h>
#include <itkAffineTransform.h>
#include <itkCompositeTransform.h>
#include <itkShrinkImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkImageMaskSpatialObject.h>
#include "MattesMutualInformation.h"
#include "RegularStepGradientDescentOptimizer.h"
#include "ConfigManager.h"

/**
 * @brief 图像配准主类
 * 
 * 整合互信息度量和梯度下降优化器,实现3D图像配准
 * 支持多分辨率金字塔策略
 * 
 * 功能特性:
 * - 支持刚体变换(6参数)和仿射变换(12参数)
 * - 支持加载初始变换(.h5文件)
 * - 支持从JSON配置文件加载参数
 * 
 * 扩展性说明:
 * - 可以修改InitializeTransform()使用不同的初始化策略
 * - 可以在Update()中添加更多预处理步骤
 * - 可以替换为不同的度量和优化器
 */
class ImageRegistration
{
public:
    using ImageType = itk::Image<float, 3>;
    using MaskImageType = itk::Image<unsigned char, 3>;
    using MaskSpatialObjectType = itk::ImageMaskSpatialObject<3>;
    using RigidTransformType = itk::Euler3DTransform<double>;
    using AffineTransformType = itk::AffineTransform<double, 3>;
    using CompositeTransformType = itk::CompositeTransform<double, 3>;
    
    // 为了兼容性保留TransformType别名 (默认刚体)
    using TransformType = RigidTransformType;
    using TransformPointer = TransformType::Pointer;
    
    using ParametersType = std::vector<double>;
    using ObserverCallbackType = std::function<void(int, double, double)>;
    using LevelObserverCallbackType = std::function<void(int, unsigned int, double)>;

    ImageRegistration();
    ~ImageRegistration();

    // =========== 图像设置 ===========
    void SetFixedImagePath(const std::string& path);
    void SetMovingImagePath(const std::string& path);
    void SetFixedImage(ImageType::Pointer fixedImage);
    void SetMovingImage(ImageType::Pointer movingImage);
    
    // =========== 掩膜设置 (用于局部配准) ===========
    // 从文件加载固定图像掩膜 (支持 .nrrd/.nii.gz 等格式)
    // 掩膜区域值>0的区域将参与配准计算,值=0的区域被忽略
    bool LoadFixedMask(const std::string& maskFilePath);
    bool HasFixedMask() const { return m_FixedImageMask.IsNotNull(); }

    // =========== 变换类型设置 ===========
    void SetTransformType(ConfigManager::TransformType type);
    ConfigManager::TransformType GetTransformType() const { return m_TransformType; }
    
    // =========== 初始变换加载 ===========
    // 从.h5文件加载初始变换(粗配准结果)
    bool LoadInitialTransform(const std::string& h5FilePath);
    void SetUseInitialTransform(bool use) { m_UseInitialTransform = use; }
    
    // =========== 从配置文件加载 ===========
    void LoadFromConfig(const ConfigManager::RegistrationConfig& config);

    // =========== 配准参数设置 ===========
    void SetNumberOfHistogramBins(unsigned int bins) { m_NumberOfHistogramBins = bins; }
    void SetNumberOfSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; }
    void SetNumberOfSpatialSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; }
    void SetLearningRate(const std::vector<double>& rates) { m_LearningRate = rates; }
    void SetMinimumStepLength(double stepLength) { m_MinimumStepLength = stepLength; }
    void SetNumberOfIterations(const std::vector<unsigned int>& iterations) { m_NumberOfIterations = iterations; }
    void SetRelaxationFactor(double factor) { m_RelaxationFactor = factor; }
    void SetGradientMagnitudeTolerance(double tol) { m_GradientMagnitudeTolerance = tol; }
    
    // =========== 多分辨率设置 ===========
    void SetNumberOfLevels(unsigned int levels) { m_NumberOfLevels = levels; }
    void SetShrinkFactors(const std::vector<unsigned int>& factors) { m_ShrinkFactors = factors; }
    void SetShrinkFactorsPerLevel(const std::vector<unsigned int>& factors) { m_ShrinkFactors = factors; }
    void SetSmoothingSigmas(const std::vector<double>& sigmas) { m_SmoothingSigmas = sigmas; }
    void SetSmoothingSigmasPerLevel(const std::vector<double>& sigmas) { m_SmoothingSigmas = sigmas; }
    
    // =========== 观察者回调 ===========
    void SetIterationObserver(ObserverCallbackType callback) { m_IterationObserver = callback; }
    void SetLevelObserver(LevelObserverCallbackType callback) { m_LevelObserver = callback; }
    
    // =========== 其他设置 ===========
    void SetRandomSeed(unsigned int seed) { m_RandomSeed = seed; }
    void SetUseStratifiedSampling(bool use) { m_UseStratifiedSampling = use; }
    void SetSamplingPercentage(double percent) { m_SamplingPercentage = percent; }
    void SetVerbose(bool v) { m_Verbose = v; }
    bool GetVerbose() const { return m_Verbose; }
    
    // =========== 参数获取 (用于级联配准) ===========
    ImageType::Pointer GetFixedImage() const { return m_FixedImage; }
    ImageType::Pointer GetMovingImage() const { return m_MovingImage; }
    MaskSpatialObjectType::Pointer GetFixedImageMask() const { return m_FixedImageMask; }
    unsigned int GetNumberOfHistogramBins() const { return m_NumberOfHistogramBins; }
    double GetSamplingPercentage() const { return m_SamplingPercentage; }
    std::vector<double> GetLearningRate() const { return m_LearningRate; }
    double GetMinimumStepLength() const { return m_MinimumStepLength; }
    std::vector<unsigned int> GetNumberOfIterations() const { return m_NumberOfIterations; }
    double GetRelaxationFactor() const { return m_RelaxationFactor; }
    double GetGradientMagnitudeTolerance() const { return m_GradientMagnitudeTolerance; }
    unsigned int GetNumberOfLevels() const { return m_NumberOfLevels; }
    std::vector<unsigned int> GetShrinkFactors() const { return m_ShrinkFactors; }
    std::vector<double> GetSmoothingSigmas() const { return m_SmoothingSigmas; }
    unsigned int GetRandomSeed() const { return m_RandomSeed; }
    void SetFixedImageMask(MaskSpatialObjectType::Pointer mask) { m_FixedImageMask = mask; }

    // =========== 执行配准 ===========
    void Update();

    // =========== 获取结果 ===========
    // 获取最终的复合变换 (包含初始变换 + 优化后的变换)
    CompositeTransformType::Pointer GetFinalCompositeTransform() const;
    
    // 获取刚体变换结果 (仅当使用刚体变换时有效)
    RigidTransformType::Pointer GetRigidTransform() const { return m_RigidTransform; }
    
    // 获取仿射变换结果 (仅当使用仿射变换时有效)
    AffineTransformType::Pointer GetAffineTransform() const { return m_AffineTransform; }
    
    // 为了兼容性保留的接口
    TransformType::Pointer GetTransform() const { return m_RigidTransform; }
    TransformType::Pointer GetFinalTransform() const { return m_RigidTransform; }

    // 获取最终度量值
    double GetFinalMetricValue() const { return m_FinalMetricValue; }
    
    // 获取配准耗时
    double GetElapsedTime() const { return m_ElapsedTime; }
    
    // 获取优化参数数量
    unsigned int GetNumberOfParameters() const;

private:
    // =========== 输入图像 ===========
    ImageType::Pointer m_FixedImage;
    ImageType::Pointer m_MovingImage;
    
    // =========== 掩膜 (用于局部配准) ===========
    MaskSpatialObjectType::Pointer m_FixedImageMask;
    unsigned long m_MaskVoxelCount;  // 掩膜内体素数 (用于正确显示采样信息)

    // =========== 变换 ===========
    ConfigManager::TransformType m_TransformType;
    RigidTransformType::Pointer m_RigidTransform;
    AffineTransformType::Pointer m_AffineTransform;
    CompositeTransformType::Pointer m_InitialTransform;
    bool m_UseInitialTransform;

    // =========== 度量和优化器 ===========
    std::unique_ptr<MattesMutualInformation> m_Metric;
    std::unique_ptr<RegularStepGradientDescentOptimizer> m_Optimizer;

    // =========== 配准参数 ===========
    unsigned int m_NumberOfHistogramBins;
    unsigned int m_NumberOfSpatialSamples;
    double m_SamplingPercentage; // 比例形式: 0.1 = 10%
    std::vector<double> m_LearningRate;  // 支持分层学习率
    double m_MinimumStepLength;
    std::vector<unsigned int> m_NumberOfIterations;  // 支持分层迭代次数
    double m_RelaxationFactor;
    double m_GradientMagnitudeTolerance;
    
    // =========== 多分辨率参数 ===========
    unsigned int m_NumberOfLevels;
    std::vector<unsigned int> m_ShrinkFactors;
    std::vector<double> m_SmoothingSigmas;
    unsigned int m_RandomSeed;
    bool m_UseStratifiedSampling;
    
    // =========== 观察者回调 ===========
    ObserverCallbackType m_IterationObserver;
    LevelObserverCallbackType m_LevelObserver;

    // =========== 结果 ===========
    double m_FinalMetricValue;
    double m_ElapsedTime;
    bool m_Verbose;

    // =========== 内部方法 ===========
    void InitializeTransform();
    void InitializeRigidTransform();
    void InitializeAffineTransform();
    
    void RunSingleLevel(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level);
    void RunSingleLevelRigid(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level);
    void RunSingleLevelAffine(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level);
    
    ImageType::Pointer ShrinkImage(ImageType::Pointer image, unsigned int factor);
    ImageType::Pointer SmoothImage(ImageType::Pointer image, double sigma);
    ImageType::Pointer WinsorizeImage(ImageType::Pointer image, double lowerQuantile = 0.005, double upperQuantile = 0.995);
    
    // 计算图像几何中心
    void ComputeGeometricCenter(ImageType::Pointer image, ImageType::PointType& center);
    
    // 自动估算参数尺度
    std::vector<double> EstimateParameterScales();
    std::vector<double> EstimateRigidParameterScales();
    std::vector<double> EstimateAffineParameterScales();
    double ComputePhysicalRadius(ImageType::Pointer image);
};

#endif // IMAGEREGISTRATION_H

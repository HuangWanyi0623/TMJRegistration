#ifndef IMAGEREGISTRATION_H
#define IMAGEREGISTRATION_H

#include <vector>
#include <memory>
#include <functional>
#include <string>
#include <itkImage.h>
#include <itkEuler3DTransform.h>
#include <itkShrinkImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include "MattesMutualInformation.h"
#include "RegularStepGradientDescentOptimizer.h"

/**
 * @brief 图像配准主类
 * 
 * 整合互信息度量和梯度下降优化器,实现3D图像配准
 * 支持多分辨率金字塔策略
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
    using TransformType = itk::Euler3DTransform<double>;
    using TransformPointer = TransformType::Pointer;
    using ParametersType = std::vector<double>;
    using ObserverCallbackType = std::function<void(int, double, double)>;
    using LevelObserverCallbackType = std::function<void(int, unsigned int, double)>;

    ImageRegistration();
    ~ImageRegistration();

    // 设置图像 (从路径加载)
    void SetFixedImagePath(const std::string& path);
    void SetMovingImagePath(const std::string& path);
    
    // 设置图像 (直接设置指针)
    void SetFixedImage(ImageType::Pointer fixedImage);
    void SetMovingImage(ImageType::Pointer movingImage);

    // 设置配准参数
    void SetNumberOfHistogramBins(unsigned int bins) { m_NumberOfHistogramBins = bins; }
    void SetNumberOfSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; }
    void SetNumberOfSpatialSamples(unsigned int samples) { m_NumberOfSpatialSamples = samples; }
    void SetLearningRate(double rate) { m_LearningRate = rate; }
    void SetMinimumStepLength(double stepLength) { m_MinimumStepLength = stepLength; }
    void SetNumberOfIterations(unsigned int iterations) { m_NumberOfIterations = iterations; }
    void SetRelaxationFactor(double factor) { m_RelaxationFactor = factor; }
    void SetGradientMagnitudeTolerance(double tol) { m_GradientMagnitudeTolerance = tol; }
    
    // 多分辨率设置
    void SetNumberOfLevels(unsigned int levels) { m_NumberOfLevels = levels; }
    void SetShrinkFactors(const std::vector<unsigned int>& factors) { m_ShrinkFactors = factors; }
    void SetShrinkFactorsPerLevel(const std::vector<unsigned int>& factors) { m_ShrinkFactors = factors; }
    void SetSmoothingSigmas(const std::vector<double>& sigmas) { m_SmoothingSigmas = sigmas; }
    void SetSmoothingSigmasPerLevel(const std::vector<double>& sigmas) { m_SmoothingSigmas = sigmas; }
    
    // 观察者回调
    void SetIterationObserver(ObserverCallbackType callback) { m_IterationObserver = callback; }
    void SetLevelObserver(LevelObserverCallbackType callback) { m_LevelObserver = callback; }
    
    // 随机种子(确保可重复性)
    void SetRandomSeed(unsigned int seed) { m_RandomSeed = seed; }

    // 执行配准
    void Update();

    // 获取变换
    TransformType::Pointer GetTransform() const { return m_Transform; }
    TransformType::Pointer GetFinalTransform() const { return m_Transform; }

    // 获取最终度量值
    double GetFinalMetricValue() const { return m_FinalMetricValue; }
    
    // 获取配准耗时
    double GetElapsedTime() const { return m_ElapsedTime; }

private:
    // 输入图像
    ImageType::Pointer m_FixedImage;
    ImageType::Pointer m_MovingImage;

    // 变换
    TransformType::Pointer m_Transform;

    // 度量和优化器
    std::unique_ptr<MattesMutualInformation> m_Metric;
    std::unique_ptr<RegularStepGradientDescentOptimizer> m_Optimizer;

    // 配准参数
    unsigned int m_NumberOfHistogramBins;
    unsigned int m_NumberOfSpatialSamples;
    double m_LearningRate;
    double m_MinimumStepLength;
    unsigned int m_NumberOfIterations;
    double m_RelaxationFactor;
    double m_GradientMagnitudeTolerance;
    
    // 多分辨率参数
    unsigned int m_NumberOfLevels;
    std::vector<unsigned int> m_ShrinkFactors;
    std::vector<double> m_SmoothingSigmas;
    unsigned int m_RandomSeed;
    
    // 观察者回调
    ObserverCallbackType m_IterationObserver;
    LevelObserverCallbackType m_LevelObserver;

    // 结果
    double m_FinalMetricValue;
    double m_ElapsedTime;

    // 内部方法
    void InitializeTransform();
    void RunSingleLevel(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, unsigned int level);
    ImageType::Pointer ShrinkImage(ImageType::Pointer image, unsigned int factor);
    ImageType::Pointer SmoothImage(ImageType::Pointer image, double sigma);
    
    // 计算图像几何中心
    void ComputeGeometricCenter(ImageType::Pointer image, ImageType::PointType& center);
};

#endif // IMAGEREGISTRATION_H

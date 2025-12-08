#ifndef REGULAR_STEP_GRADIENT_DESCENT_OPTIMIZER_H
#define REGULAR_STEP_GRADIENT_DESCENT_OPTIMIZER_H

#include <vector>
#include <functional>

/**
 * @brief 规则步长梯度下降优化器
 * 
 * 基于ITK的RegularStepGradientDescentOptimizerv4实现
 * 关键特性:
 * - 参数尺度(Scales): 处理旋转(rad)和平移(mm)的量级差异
 * - 回退机制: 代价上升时回退参数并减小步长
 * - 与ITK行为一致的梯度处理
 * 
 * 扩展性说明:
 * - 可以修改AdvanceOneStep()改变更新策略
 * - 可以添加动量项实现加速收敛
 * - 可以实现自适应学习率策略
 */
class RegularStepGradientDescentOptimizer
{
public:
    using ParametersType = std::vector<double>;
    using CostFunctionType = std::function<double()>;
    using GradientFunctionType = std::function<void(ParametersType&)>;
    using UpdateParametersType = std::function<void(const ParametersType&)>;
    using GetParametersType = std::function<ParametersType()>;
    using SetParametersType = std::function<void(const ParametersType&)>;
    using ObserverType = std::function<void(unsigned int, double, double)>;

    RegularStepGradientDescentOptimizer();
    ~RegularStepGradientDescentOptimizer();

    // 设置优化参数
    void SetLearningRate(double rate) { m_LearningRate = rate; m_CurrentStepLength = rate; }
    void SetMinimumStepLength(double stepLength) { m_MinimumStepLength = stepLength; }
    void SetNumberOfIterations(unsigned int iterations) { m_NumberOfIterations = iterations; }
    void SetRelaxationFactor(double factor) { m_RelaxationFactor = factor; }
    void SetGradientMagnitudeTolerance(double tolerance) { m_GradientMagnitudeTolerance = tolerance; }
    void SetReturnBestParametersAndValue(bool flag) { m_ReturnBestParameters = flag; }
    
    // 设置参数尺度 - 关键! 用于平衡旋转和平移参数的更新幅度
    void SetScales(const ParametersType& scales) { m_Scales = scales; }

    // 设置代价函数和梯度函数
    void SetCostFunction(CostFunctionType costFunc) { m_CostFunction = costFunc; }
    void SetGradientFunction(GradientFunctionType gradFunc) { m_GradientFunction = gradFunc; }
    void SetUpdateParametersFunction(UpdateParametersType updateFunc) { m_UpdateParameters = updateFunc; }
    
    // 设置参数获取和设置函数 (用于回退机制)
    void SetGetParametersFunction(GetParametersType getFunc) { m_GetParameters = getFunc; }
    void SetSetParametersFunction(SetParametersType setFunc) { m_SetParameters = setFunc; }

    // 设置参数数量
    void SetNumberOfParameters(unsigned int num);

    // 执行优化
    void StartOptimization();

    // 获取结果
    double GetValue() const { return m_CurrentValue; }
    double GetBestValue() const { return m_BestValue; }
    unsigned int GetCurrentIteration() const { return m_CurrentIteration; }
    double GetLearningRate() const { return m_CurrentStepLength; }
    
    // 停止原因
    enum StopCondition { MAXIMUM_ITERATIONS, STEP_TOO_SMALL, GRADIENT_TOO_SMALL, CONVERGED };
    StopCondition GetStopCondition() const { return m_StopCondition; }

    // 设置观察者回调(用于输出优化过程)
    void SetObserver(ObserverType observer) { m_Observer = observer; }

    // 设置观察者每隔多少迭代调用一次 (默认 10), 如果 verbose 则每次迭代调用
    void SetObserverIterationInterval(unsigned int interval) { m_ObserverIterationInterval = interval; }
    unsigned int GetObserverIterationInterval() const { return m_ObserverIterationInterval; }

    // 调试日志
    void SetVerbose(bool verbose) { m_Verbose = verbose; }
    bool GetVerbose() const { return m_Verbose; }

    // 设置参数更新的最大值以避免过大的跳跃(单位根据参数：弧度或mm)
    void SetMaxParameterUpdate(const ParametersType& maxUpdate) { m_MaxParameterUpdate = maxUpdate; }
    const ParametersType& GetMaxParameterUpdate() const { return m_MaxParameterUpdate; }

private:
    // 优化参数
    double m_LearningRate;
    double m_MinimumStepLength;
    unsigned int m_NumberOfIterations;
    double m_RelaxationFactor;
    double m_GradientMagnitudeTolerance;
    bool m_ReturnBestParameters;
    unsigned int m_NumberOfParameters;
    ParametersType m_Scales;  // 参数尺度

    // 当前状态
    double m_CurrentValue;
    double m_BestValue;
    unsigned int m_CurrentIteration;
    double m_CurrentStepLength;
    double m_PreviousValue;
    ParametersType m_CurrentGradient;
    ParametersType m_PreviousParameters;  // 用于回退
    ParametersType m_BestParameters;      // 最佳参数
    ParametersType m_MaxParameterUpdate; // 每个参数最大更新值
    StopCondition m_StopCondition;

    // 代价函数和梯度函数
    CostFunctionType m_CostFunction;
    GradientFunctionType m_GradientFunction;
    UpdateParametersType m_UpdateParameters;
    GetParametersType m_GetParameters;
    SetParametersType m_SetParameters;

    // 观察者
    ObserverType m_Observer;
    unsigned int m_ObserverIterationInterval;
    bool m_Verbose;

    // 内部方法
    void AdvanceOneStep();
    double ComputeScaledGradientMagnitude(const ParametersType& gradient);
};

#endif // REGULAR_STEP_GRADIENT_DESCENT_OPTIMIZER_H

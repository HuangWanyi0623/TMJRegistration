#ifndef REGULAR_STEP_GRADIENT_DESCENT_OPTIMIZER_H
#define REGULAR_STEP_GRADIENT_DESCENT_OPTIMIZER_H

#include <vector>
#include <functional>

/**
 * @brief 规则步长梯度下降优化器
 * 
 * 基于ITK的RegularStepGradientDescentOptimizerv4实现
 */
class RegularStepGradientDescentOptimizer
{
public:
    using ParametersType = std::vector<double>;
    using CostFunctionType = std::function<double(const ParametersType&)>;
    using GradientFunctionType = std::function<void(const ParametersType&, ParametersType&)>;

    RegularStepGradientDescentOptimizer();
    ~RegularStepGradientDescentOptimizer();

    // 设置优化参数
    void SetMaximumStepLength(double stepLength) { m_MaximumStepLength = stepLength; }
    void SetMinimumStepLength(double stepLength) { m_MinimumStepLength = stepLength; }
    void SetNumberOfIterations(unsigned int iterations) { m_NumberOfIterations = iterations; }
    void SetRelaxationFactor(double factor) { m_RelaxationFactor = factor; }
    void SetGradientMagnitudeTolerance(double tolerance) { m_GradientMagnitudeTolerance = tolerance; }

    // 设置初始位置
    void SetInitialPosition(const ParametersType& initialPosition) { m_InitialPosition = initialPosition; }

    // 设置代价函数和梯度函数
    void SetCostFunction(CostFunctionType costFunc) { m_CostFunction = costFunc; }
    void SetGradientFunction(GradientFunctionType gradFunc) { m_GradientFunction = gradFunc; }

    // 执行优化
    void StartOptimization();

    // 获取结果
    const ParametersType& GetCurrentPosition() const { return m_CurrentPosition; }
    double GetCurrentValue() const { return m_CurrentValue; }
    unsigned int GetCurrentIteration() const { return m_CurrentIteration; }

    // 设置观察者回调(用于输出优化过程)
    void SetObserver(std::function<void(unsigned int, double, const ParametersType&)> observer) {
        m_Observer = observer;
    }

private:
    // 优化参数
    double m_MaximumStepLength;
    double m_MinimumStepLength;
    unsigned int m_NumberOfIterations;
    double m_RelaxationFactor;
    double m_GradientMagnitudeTolerance;

    // 当前状态
    ParametersType m_InitialPosition;
    ParametersType m_CurrentPosition;
    double m_CurrentValue;
    unsigned int m_CurrentIteration;
    double m_CurrentStepLength;

    // 代价函数和梯度函数
    CostFunctionType m_CostFunction;
    GradientFunctionType m_GradientFunction;

    // 观察者
    std::function<void(unsigned int, double, const ParametersType&)> m_Observer;

    // 内部方法
    void AdvanceOneStep();
    double ComputeGradientMagnitude(const ParametersType& gradient);
};

#endif // REGULAR_STEP_GRADIENT_DESCENT_OPTIMIZER_H

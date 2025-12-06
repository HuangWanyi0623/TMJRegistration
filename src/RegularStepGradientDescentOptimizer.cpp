#include "RegularStepGradientDescentOptimizer.h"
#include <cmath>
#include <iostream>
#include <iomanip>

RegularStepGradientDescentOptimizer::RegularStepGradientDescentOptimizer()
    : m_MaximumStepLength(1.0)
    , m_MinimumStepLength(0.001)
    , m_NumberOfIterations(100)
    , m_RelaxationFactor(0.5)
    , m_GradientMagnitudeTolerance(1e-4)
    , m_CurrentValue(0.0)
    , m_CurrentIteration(0)
    , m_CurrentStepLength(1.0)
{
}

RegularStepGradientDescentOptimizer::~RegularStepGradientDescentOptimizer()
{
}

void RegularStepGradientDescentOptimizer::StartOptimization()
{
    if (!m_CostFunction || !m_GradientFunction)
    {
        throw std::runtime_error("Cost function or gradient function not set");
    }

    if (m_InitialPosition.empty())
    {
        throw std::runtime_error("Initial position not set");
    }

    // 初始化
    m_CurrentPosition = m_InitialPosition;
    m_CurrentStepLength = m_MaximumStepLength;
    m_CurrentIteration = 0;

    std::cout << "\n=== Starting Optimization ===" << std::endl;
    std::cout << "Maximum iterations: " << m_NumberOfIterations << std::endl;
    std::cout << "Initial step length: " << m_MaximumStepLength << std::endl;
    std::cout << "Minimum step length: " << m_MinimumStepLength << std::endl;
    std::cout << "Relaxation factor: " << m_RelaxationFactor << std::endl;
    std::cout << "Number of parameters: " << m_CurrentPosition.size() << std::endl;

    // 开始迭代
    for (m_CurrentIteration = 0; m_CurrentIteration < m_NumberOfIterations; ++m_CurrentIteration)
    {
        // 执行一步优化
        AdvanceOneStep();

        // 调用观察者
        if (m_Observer)
        {
            m_Observer(m_CurrentIteration, m_CurrentValue, m_CurrentPosition);
        }

        // 检查收敛条件
        if (m_CurrentStepLength < m_MinimumStepLength)
        {
            std::cout << "\nOptimization converged: step length below minimum" << std::endl;
            break;
        }
    }

    std::cout << "\n=== Optimization Completed ===" << std::endl;
    std::cout << "Final iteration: " << m_CurrentIteration << std::endl;
    std::cout << "Final metric value: " << m_CurrentValue << std::endl;
    std::cout << "Final step length: " << m_CurrentStepLength << std::endl;
}

void RegularStepGradientDescentOptimizer::AdvanceOneStep()
{
    // 计算当前梯度
    ParametersType gradient(m_CurrentPosition.size());
    m_GradientFunction(m_CurrentPosition, gradient);

    // 计算梯度幅值
    double gradientMagnitude = ComputeGradientMagnitude(gradient);

    if (m_CurrentIteration % 10 == 0)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Iter " << std::setw(4) << m_CurrentIteration 
                  << ": Value = " << std::setw(10) << m_CurrentValue
                  << ", GradMag = " << std::setw(10) << gradientMagnitude
                  << ", Step = " << std::setw(8) << m_CurrentStepLength << std::endl;
    }

    // 检查梯度是否过小
    if (gradientMagnitude < m_GradientMagnitudeTolerance)
    {
        std::cout << "\nGradient magnitude below tolerance" << std::endl;
        return;
    }

    // 归一化梯度
    for (size_t i = 0; i < gradient.size(); ++i)
    {
        gradient[i] /= gradientMagnitude;
    }

    // 更新参数
    ParametersType newPosition = m_CurrentPosition;
    for (size_t i = 0; i < newPosition.size(); ++i)
    {
        newPosition[i] -= m_CurrentStepLength * gradient[i];
    }

    // 计算新位置的代价
    double newValue = m_CostFunction(newPosition);

    // 如果代价下降,接受新位置
    if (newValue < m_CurrentValue || m_CurrentIteration == 0)
    {
        m_CurrentPosition = newPosition;
        m_CurrentValue = newValue;
    }
    else
    {
        // 如果代价上升,减小步长
        m_CurrentStepLength *= m_RelaxationFactor;
    }
}

double RegularStepGradientDescentOptimizer::ComputeGradientMagnitude(const ParametersType& gradient)
{
    double magnitude = 0.0;
    for (double g : gradient)
    {
        magnitude += g * g;
    }
    return std::sqrt(magnitude);
}

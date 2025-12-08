#include "RegularStepGradientDescentOptimizer.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

RegularStepGradientDescentOptimizer::RegularStepGradientDescentOptimizer()
    : m_LearningRate(1.0)
    , m_MinimumStepLength(0.0001)
    , m_NumberOfIterations(300)
    , m_RelaxationFactor(0.5)
    , m_GradientMagnitudeTolerance(1e-4)
    , m_ReturnBestParameters(true)
    , m_NumberOfParameters(6)
    , m_CurrentValue(0.0)
    , m_BestValue(std::numeric_limits<double>::max())
    , m_CurrentIteration(0)
    , m_CurrentStepLength(1.0)
    , m_PreviousValue(std::numeric_limits<double>::max())
    , m_StopCondition(MAXIMUM_ITERATIONS)
    , m_Verbose(false)
    , m_ObserverIterationInterval(10)
    , m_MaxParameterUpdate(m_NumberOfParameters, 1.0) // default max update
{
    // 默认尺度 - 全为1
    m_Scales.resize(6, 1.0);
    // Default max parameter updates: rotations ~0.2 rad, translations ~20mm
    m_MaxParameterUpdate = ParametersType{0.2, 0.2, 0.2, 20.0, 20.0, 20.0};
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
    
    if (!m_GetParameters || !m_SetParameters)
    {
        throw std::runtime_error("Get/Set parameters functions not set");
    }

    // 初始化
    m_CurrentStepLength = m_LearningRate;
    m_CurrentIteration = 0;
    m_BestValue = std::numeric_limits<double>::max();
    m_PreviousValue = std::numeric_limits<double>::max();
    m_CurrentGradient.resize(m_NumberOfParameters, 0.0);
    
    // 确保尺度正确
    if (m_Scales.size() != m_NumberOfParameters)
    {
        m_Scales.resize(m_NumberOfParameters, 1.0);
    }

    // 获取初始参数
    m_PreviousParameters = m_GetParameters();
    m_BestParameters = m_PreviousParameters;
    
    // 计算初始值
    m_CurrentValue = m_CostFunction();
    m_BestValue = m_CurrentValue;

    if (m_Verbose)
    {
        std::cout << "[Optimizer Debug] Initial parameters: ";
        for (double p : m_PreviousParameters) std::cout << std::fixed << std::setprecision(6) << p << " ";
        std::cout << "\n";
        std::cout << "[Optimizer Debug] Initial value: " << m_CurrentValue << std::endl;
    }

    // 开始迭代
    for (m_CurrentIteration = 0; m_CurrentIteration < m_NumberOfIterations; ++m_CurrentIteration)
    {
        // 调用观察者
        if (m_Observer && (m_Verbose || (m_CurrentIteration % m_ObserverIterationInterval == 0)))
        {
            m_Observer(m_CurrentIteration, m_CurrentValue, m_CurrentStepLength);
        }

        // 执行一步优化
        AdvanceOneStep();

        // 检查收敛条件
        if (m_CurrentStepLength < m_MinimumStepLength)
        {
            m_StopCondition = STEP_TOO_SMALL;
            break;
        }
    }

    if (m_CurrentIteration >= m_NumberOfIterations)
    {
        m_StopCondition = MAXIMUM_ITERATIONS;
    }
    
    // 恢复最佳参数
    if (m_ReturnBestParameters && !m_BestParameters.empty())
    {
        m_SetParameters(m_BestParameters);
        m_CurrentValue = m_BestValue;
    }
}

void RegularStepGradientDescentOptimizer::AdvanceOneStep()
{
    // 保存当前参数(用于可能的回退)
    m_PreviousParameters = m_GetParameters();
    m_PreviousValue = m_CurrentValue;
    
    // 计算当前梯度
    m_GradientFunction(m_CurrentGradient);

    // 计算考虑尺度的梯度幅值
    double gradientMagnitude = ComputeScaledGradientMagnitude(m_CurrentGradient);

    // 检查梯度是否过小
    if (gradientMagnitude < m_GradientMagnitudeTolerance)
    {
        m_StopCondition = GRADIENT_TOO_SMALL;
        m_CurrentStepLength = 0.0; // 触发停止
        return;
    }

    // 计算缩放后的归一化梯度方向并应用步长
    // ITK的方式: direction = gradient / (scales^2 * magnitude)
    // 这样每个参数的更新量 = stepLength * gradient[i] / (scales[i]^2 * magnitude)
    std::vector<double> newParameters(m_NumberOfParameters);
    auto currentParams = m_PreviousParameters;
    
    for (unsigned int i = 0; i < m_NumberOfParameters; ++i)
    {
        double scaleFactor = m_Scales[i] * m_Scales[i];
        double direction = m_CurrentGradient[i] / (scaleFactor * gradientMagnitude);
        // 负方向,因为我们要最小化
        newParameters[i] = currentParams[i] - m_CurrentStepLength * direction;
    }

    // clamp parameter updates to avoid extreme jumps
    for (unsigned int i = 0; i < m_NumberOfParameters && i < m_MaxParameterUpdate.size(); ++i)
    {
        double delta = newParameters[i] - currentParams[i];
        if (delta > m_MaxParameterUpdate[i]) delta = m_MaxParameterUpdate[i];
        if (delta < -m_MaxParameterUpdate[i]) delta = -m_MaxParameterUpdate[i];
        newParameters[i] = currentParams[i] + delta;
    }

    // 打印调试信息
    if (m_Verbose)
    {
        std::cout << "[Optimizer Debug] Iter=" << m_CurrentIteration << " gradMag=" << gradientMagnitude
                  << " step=" << m_CurrentStepLength << "\n";
        std::cout << "  gradient: ";
        for (unsigned int i = 0; i < std::min<unsigned int>(m_CurrentGradient.size(), 6); ++i)
            std::cout << std::fixed << std::setprecision(6) << m_CurrentGradient[i] << " ";
        std::cout << "\n";
        std::cout << "  param update: ";
        for (unsigned int i = 0; i < std::min<unsigned int>(newParameters.size(), 6); ++i)
            std::cout << std::fixed << std::setprecision(6) << newParameters[i] - currentParams[i] << " ";
        std::cout << "\n";
    }

    // 应用新参数
    m_SetParameters(newParameters);

    // 计算新的代价值
    m_CurrentValue = m_CostFunction();

    // 检查是否改进
    if (m_CurrentValue < m_PreviousValue)
    {
        // 代价下降,接受这步
        if (m_CurrentValue < m_BestValue)
        {
            m_BestValue = m_CurrentValue;
            m_BestParameters = newParameters;
        }
    }
    else
    {
        // 代价上升或相等,回退参数并减小步长
        m_SetParameters(m_PreviousParameters);
        m_CurrentValue = m_PreviousValue;
        m_CurrentStepLength *= m_RelaxationFactor;
    }
}

double RegularStepGradientDescentOptimizer::ComputeScaledGradientMagnitude(const ParametersType& gradient)
{
    // ITK计算方式: sqrt(sum( (gradient[i] / scale[i])^2 ))
    double magnitude = 0.0;
    for (unsigned int i = 0; i < gradient.size() && i < m_Scales.size(); ++i)
    {
        double scaledGrad = gradient[i] / m_Scales[i];
        magnitude += scaledGrad * scaledGrad;
    }
    return std::sqrt(magnitude);
}

void RegularStepGradientDescentOptimizer::SetNumberOfParameters(unsigned int num)
{
    m_NumberOfParameters = num;
    // resize gradient and max updates if needed
    m_CurrentGradient.resize(m_NumberOfParameters);
    if (m_MaxParameterUpdate.size() != m_NumberOfParameters)
    {
        // keep existing values or set sensible defaults
        m_MaxParameterUpdate.resize(m_NumberOfParameters, 1.0);
        if (m_NumberOfParameters >= 6)
        {
            // default for rigid-like: rotation small, translation bigger
            for (unsigned int i = 0; i < m_NumberOfParameters; ++i)
            {
                if (i < 3) m_MaxParameterUpdate[i] = 0.2; else m_MaxParameterUpdate[i] = 20.0;
            }
        }
    }
}

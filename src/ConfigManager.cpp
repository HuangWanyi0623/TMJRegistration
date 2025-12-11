#include "ConfigManager.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>

// ============================================================================
// 构造函数和析构函数
// ============================================================================

ConfigManager::ConfigManager()
{
    // 使用默认配置
}

ConfigManager::~ConfigManager()
{
}

// ============================================================================
// 变换类型转换
// ============================================================================

std::string ConfigManager::TransformTypeToString(TransformType type)
{
    switch (type)
    {
        case TransformType::Rigid: return "Rigid";
        case TransformType::Affine: return "Affine";
        case TransformType::RigidThenAffine: return "RigidThenAffine";
        default: return "Rigid";
    }
}

ConfigManager::TransformType ConfigManager::StringToTransformType(const std::string& str)
{
    std::string lower = str;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    if (lower == "affine") return TransformType::Affine;
    if (lower == "rigidthenaffine" || lower == "rigid+affine" || lower == "rigidaffine") 
        return TransformType::RigidThenAffine;
    return TransformType::Rigid;  // 默认刚体
}

void ConfigManager::SetTransformType(const std::string& typeStr)
{
    m_Config.transformType = StringToTransformType(typeStr);
}

// ============================================================================
// 字符串处理辅助函数
// ============================================================================

std::string ConfigManager::Trim(const std::string& str) const
{
    size_t start = str.find_first_not_of(" \t\n\r\"");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\n\r\"");
    return str.substr(start, end - start + 1);
}

// ============================================================================
// 简单JSON解析
// ============================================================================

std::string ConfigManager::ExtractValue(const std::string& content, const std::string& key) const
{
    // 查找 "key": value 或 "key": "value"
    std::string searchKey = "\"" + key + "\"";
    size_t pos = content.find(searchKey);
    if (pos == std::string::npos) return "";
    
    pos = content.find(':', pos);
    if (pos == std::string::npos) return "";
    
    pos++;
    
    // 跳过空白
    while (pos < content.size() && std::isspace(content[pos])) pos++;
    
    if (pos >= content.size()) return "";
    
    // 检查是否是数组
    if (content[pos] == '[')
    {
        return "";  // 使用ExtractArray处理数组
    }
    
    // 查找值的结束位置 (逗号, 右花括号, 或换行)
    size_t endPos = content.find_first_of(",}\n", pos);
    if (endPos == std::string::npos) endPos = content.size();
    
    return Trim(content.substr(pos, endPos - pos));
}

std::vector<std::string> ConfigManager::ExtractArray(const std::string& content, const std::string& key) const
{
    std::vector<std::string> result;
    
    std::string searchKey = "\"" + key + "\"";
    size_t pos = content.find(searchKey);
    if (pos == std::string::npos) return result;
    
    pos = content.find('[', pos);
    if (pos == std::string::npos) return result;
    
    size_t endPos = content.find(']', pos);
    if (endPos == std::string::npos) return result;
    
    std::string arrayContent = content.substr(pos + 1, endPos - pos - 1);
    
    // 分割数组元素
    std::stringstream ss(arrayContent);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        std::string trimmed = Trim(item);
        if (!trimmed.empty())
        {
            result.push_back(trimmed);
        }
    }
    
    return result;
}

bool ConfigManager::ParseJsonFile(const std::string& content)
{
    try
    {
        // 解析变换类型
        std::string transformType = ExtractValue(content, "transformType");
        if (!transformType.empty())
        {
            m_Config.transformType = StringToTransformType(transformType);
        }
        
        // 解析度量参数
        std::string bins = ExtractValue(content, "numberOfHistogramBins");
        if (!bins.empty()) m_Config.numberOfHistogramBins = std::stoul(bins);
        
    std::string samples = ExtractValue(content, "numberOfSpatialSamples");
    if (!samples.empty()) m_Config.numberOfSpatialSamples = std::stoul(samples);
    std::string sampPct = ExtractValue(content, "samplingPercentage");
    if (!sampPct.empty()) m_Config.samplingPercentage = std::stod(sampPct);
        
        // 解析优化器参数 - 学习率支持单值或数组
        std::string lr = ExtractValue(content, "learningRate");
        if (!lr.empty())
        {
            // 尝试作为单个值解析
            try
            {
                double singleLR = std::stod(lr);
                m_Config.learningRate.clear();
                m_Config.learningRate.push_back(singleLR);
            }
            catch (...)
            {
                // 如果失败，尝试作为数组解析
                auto lrArray = ExtractArray(content, "learningRate");
                if (!lrArray.empty())
                {
                    m_Config.learningRate.clear();
                    for (const auto& s : lrArray)
                    {
                        m_Config.learningRate.push_back(std::stod(s));
                    }
                }
            }
        }
        else
        {
            // 尝试作为数组解析
            auto lrArray = ExtractArray(content, "learningRate");
            if (!lrArray.empty())
            {
                m_Config.learningRate.clear();
                for (const auto& s : lrArray)
                {
                    m_Config.learningRate.push_back(std::stod(s));
                }
            }
        }
        
        std::string minStep = ExtractValue(content, "minimumStepLength");
        if (!minStep.empty()) m_Config.minimumStepLength = std::stod(minStep);
        
        // 解析迭代次数 - 支持单值或数组
        std::string iter = ExtractValue(content, "numberOfIterations");
        if (!iter.empty())
        {
            // 尝试作为单个值解析
            try
            {
                unsigned int singleIter = std::stoul(iter);
                m_Config.numberOfIterations.clear();
                m_Config.numberOfIterations.push_back(singleIter);
            }
            catch (...)
            {
                // 如果失败，尝试作为数组解析
                auto iterArray = ExtractArray(content, "numberOfIterations");
                if (!iterArray.empty())
                {
                    m_Config.numberOfIterations.clear();
                    for (const auto& s : iterArray)
                    {
                        m_Config.numberOfIterations.push_back(std::stoul(s));
                    }
                }
            }
        }
        else
        {
            // 尝试作为数组解析
            auto iterArray = ExtractArray(content, "numberOfIterations");
            if (!iterArray.empty())
            {
                m_Config.numberOfIterations.clear();
                for (const auto& s : iterArray)
                {
                    m_Config.numberOfIterations.push_back(std::stoul(s));
                }
            }
        }
        
        std::string relax = ExtractValue(content, "relaxationFactor");
        if (!relax.empty()) m_Config.relaxationFactor = std::stod(relax);
        
        std::string gradTol = ExtractValue(content, "gradientMagnitudeTolerance");
        if (!gradTol.empty()) m_Config.gradientMagnitudeTolerance = std::stod(gradTol);
        
        // 解析多分辨率参数
        std::string levels = ExtractValue(content, "numberOfLevels");
        if (!levels.empty()) m_Config.numberOfLevels = std::stoul(levels);
        
        auto shrinkArray = ExtractArray(content, "shrinkFactors");
        if (!shrinkArray.empty())
        {
            m_Config.shrinkFactors.clear();
            for (const auto& s : shrinkArray)
            {
                m_Config.shrinkFactors.push_back(std::stoul(s));
            }
        }
        
        auto sigmaArray = ExtractArray(content, "smoothingSigmas");
        if (!sigmaArray.empty())
        {
            m_Config.smoothingSigmas.clear();
            for (const auto& s : sigmaArray)
            {
                m_Config.smoothingSigmas.push_back(std::stod(s));
            }
        }
        
        // 解析采样参数
        std::string stratified = ExtractValue(content, "useStratifiedSampling");
        if (!stratified.empty())
        {
            std::string lower = stratified;
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            m_Config.useStratifiedSampling = (lower == "true" || lower == "1");
        }
        
        std::string seed = ExtractValue(content, "randomSeed");
        if (!seed.empty()) m_Config.randomSeed = std::stoul(seed);
        
        return true;
    }
    catch (const std::exception& e)
    {
        std::cerr << "[Config Error] Failed to parse config: " << e.what() << std::endl;
        return false;
    }
}

// ============================================================================
// 文件I/O
// ============================================================================

bool ConfigManager::LoadFromFile(const std::string& filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open())
    {
        std::cerr << "[Config] Could not open config file: " << filePath << std::endl;
        std::cerr << "[Config] Using default configuration." << std::endl;
        return false;
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();
    file.close();
    
    bool success = ParseJsonFile(content);
    if (success)
    {
        std::cout << "[Config] Loaded configuration from: " << filePath << std::endl;
    }
    
    return success;
}

std::string ConfigManager::GenerateJson() const
{
    std::ostringstream oss;
    oss << "{\n";
    oss << "    \"_comment\": \"Registration Configuration File\",\n";
    oss << "    \n";
    oss << "    \"transformType\": \"" << TransformTypeToString(m_Config.transformType) << "\",\n";
    oss << "    \n";
    oss << "    \"_section_metric\": \"=== Metric Parameters ===\",\n";
    oss << "    \"numberOfHistogramBins\": " << m_Config.numberOfHistogramBins << ",\n";
    if (m_Config.numberOfSpatialSamples > 0)
    {
        oss << "    \"numberOfSpatialSamples\": " << m_Config.numberOfSpatialSamples << ",\n";
    }
    oss << "    \"samplingPercentage\": " << std::fixed << std::setprecision(3) << m_Config.samplingPercentage << ",\n";
    oss << "    \n";
    oss << "    \"_section_optimizer\": \"=== Optimizer Parameters ===\",\n";
    
    // 输出学习率数组
    oss << "    \"learningRate\": [";
    for (size_t i = 0; i < m_Config.learningRate.size(); ++i)
    {
        oss << std::fixed << std::setprecision(4) << m_Config.learningRate[i];
        if (i < m_Config.learningRate.size() - 1) oss << ", ";
    }
    oss << "],\n";
    
    oss << "    \"minimumStepLength\": " << std::scientific << std::setprecision(4) << m_Config.minimumStepLength << ",\n";
    
    // 输出迭代次数数组
    oss << "    \"numberOfIterations\": [";
    for (size_t i = 0; i < m_Config.numberOfIterations.size(); ++i)
    {
        oss << m_Config.numberOfIterations[i];
        if (i < m_Config.numberOfIterations.size() - 1) oss << ", ";
    }
    oss << "],\n";
    
    oss << "    \"relaxationFactor\": " << std::fixed << std::setprecision(2) << m_Config.relaxationFactor << ",\n";
    oss << "    \"gradientMagnitudeTolerance\": " << std::scientific << std::setprecision(1) << m_Config.gradientMagnitudeTolerance << ",\n";
    oss << "    \n";
    oss << "    \"_section_multiresolution\": \"=== Multi-Resolution Parameters ===\",\n";
    oss << "    \"numberOfLevels\": " << m_Config.numberOfLevels << ",\n";
    
    // shrinkFactors数组
    oss << "    \"shrinkFactors\": [";
    for (size_t i = 0; i < m_Config.shrinkFactors.size(); ++i)
    {
        oss << m_Config.shrinkFactors[i];
        if (i < m_Config.shrinkFactors.size() - 1) oss << ", ";
    }
    oss << "],\n";
    
    // smoothingSigmas数组
    oss << "    \"smoothingSigmas\": [";
    for (size_t i = 0; i < m_Config.smoothingSigmas.size(); ++i)
    {
        oss << std::fixed << std::setprecision(1) << m_Config.smoothingSigmas[i];
        if (i < m_Config.smoothingSigmas.size() - 1) oss << ", ";
    }
    oss << "],\n";
    
    oss << "    \n";
    oss << "    \"_section_sampling\": \"=== Sampling Parameters ===\",\n";
    oss << "    \"useStratifiedSampling\": " << (m_Config.useStratifiedSampling ? "true" : "false") << ",\n";
    oss << "    \"randomSeed\": " << m_Config.randomSeed << "\n";
    oss << "}\n";
    
    return oss.str();
}

bool ConfigManager::SaveToFile(const std::string& filePath) const
{
    std::ofstream file(filePath);
    if (!file.is_open())
    {
        std::cerr << "[Config Error] Could not create config file: " << filePath << std::endl;
        return false;
    }
    
    file << GenerateJson();
    file.close();
    
    std::cout << "[Config] Configuration saved to: " << filePath << std::endl;
    return true;
}

bool ConfigManager::CreateDefaultConfigFile(const std::string& filePath, TransformType type)
{
    ConfigManager config;
    config.m_Config.transformType = type;
    return config.SaveToFile(filePath);
}

// ============================================================================
// 打印配置
// ============================================================================

void ConfigManager::PrintConfig() const
{
    std::cout << "\n[Configuration]" << std::endl;
    std::cout << "  Transform Type: " << TransformTypeToString(m_Config.transformType) << std::endl;
    std::cout << "  Histogram Bins: " << m_Config.numberOfHistogramBins << std::endl;
    std::cout << "  Spatial Samples: " << m_Config.numberOfSpatialSamples << std::endl;
    std::cout << "  Sampling Percentage: " << m_Config.samplingPercentage << std::endl;
    
    // 打印学习率数组
    std::cout << "  Learning Rate: [";
    for (size_t i = 0; i < m_Config.learningRate.size(); ++i)
    {
        std::cout << m_Config.learningRate[i];
        if (i < m_Config.learningRate.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    std::cout << "  Min Step Length: " << m_Config.minimumStepLength << std::endl;
    
    // 打印迭代次数数组
    std::cout << "  Max Iterations: [";
    for (size_t i = 0; i < m_Config.numberOfIterations.size(); ++i)
    {
        std::cout << m_Config.numberOfIterations[i];
        if (i < m_Config.numberOfIterations.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    std::cout << "  Relaxation Factor: " << m_Config.relaxationFactor << std::endl;
    std::cout << "  Gradient Tolerance: " << m_Config.gradientMagnitudeTolerance << std::endl;
    std::cout << "  Multi-Resolution Levels: " << m_Config.numberOfLevels << std::endl;
    
    std::cout << "  Shrink Factors: [";
    for (size_t i = 0; i < m_Config.shrinkFactors.size(); ++i)
    {
        std::cout << m_Config.shrinkFactors[i];
        if (i < m_Config.shrinkFactors.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    std::cout << "  Smoothing Sigmas: [";
    for (size_t i = 0; i < m_Config.smoothingSigmas.size(); ++i)
    {
        std::cout << m_Config.smoothingSigmas[i];
        if (i < m_Config.smoothingSigmas.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    std::cout << "  Stratified Sampling: " << (m_Config.useStratifiedSampling ? "Yes" : "No") << std::endl;
    std::cout << "  Random Seed: " << m_Config.randomSeed << std::endl;
}

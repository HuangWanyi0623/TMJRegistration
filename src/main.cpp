#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

#include "ImageRegistration.h"
#include "ConfigManager.h"

#include <itkTransformFileWriter.h>
#include <itkCompositeTransform.h>

namespace fs = std::filesystem;

// ============================================================================
// 命令行参数结构
// ============================================================================

struct CommandLineArgs
{
    std::string fixedImagePath;
    std::string movingImagePath;
    std::string outputFolder;
    std::string configFilePath;
    std::string initialTransformPath;
    std::string transformType = "Rigid";  // 默认刚体
    bool showHelp = false;
    bool generateConfig = false;
    double samplingPercentage = -1.0;
    bool verbose = false;
};

// ============================================================================
// 辅助函数
// ============================================================================

#ifdef _WIN32
std::string WideToUtf8(const std::wstring& wstr)
{
    if (wstr.empty()) return std::string();
    int size_needed = WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), 
                                          static_cast<int>(wstr.size()), 
                                          nullptr, 0, nullptr, nullptr);
    std::string strTo(size_needed, 0);
    WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), 
                        static_cast<int>(wstr.size()), 
                        &strTo[0], size_needed, nullptr, nullptr);
    return strTo;
}
#endif

std::string GenerateTimestampFilename()
{
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::tm tm_now;
#ifdef _WIN32
    localtime_s(&tm_now, &time_t_now);
#else
    localtime_r(&time_t_now, &tm_now);
#endif
    
    std::ostringstream oss;
    oss << "registration_transform_"
        << std::put_time(&tm_now, "%Y%m%d_%H%M%S")
        << ".h5";
    return oss.str();
}

// ============================================================================
// 保存变换 (已废弃 - 现在直接在main中内联保存)
// ============================================================================

// bool SaveTransform(...) - 已移除，直接在main中实现

// ============================================================================
// 帮助信息
// ============================================================================

void PrintUsage(const char* programName)
{
    std::cout << "\n=== Mutual Information Registration Tool ===" << std::endl;
    std::cout << "Usage: " << programName << " [options] <fixed_image> <moving_image> <output_folder>\n" << std::endl;
    
    std::cout << "Required Arguments:" << std::endl;
    std::cout << "  <fixed_image>       Path to the fixed (reference) image" << std::endl;
    std::cout << "  <moving_image>      Path to the moving image to be registered" << std::endl;
    std::cout << "  <output_folder>     Folder path for output .h5 transform file\n" << std::endl;
    
    std::cout << "Options:" << std::endl;
    std::cout << "  --help, -h          Show this help message" << std::endl;
    std::cout << "  --config <file>     Load configuration from JSON file" << std::endl;
    std::cout << "  --initial <file>    Load initial transform from .h5 file (coarse registration)" << std::endl;
    std::cout << "  --transform <type>  Transform type: Rigid (default) or Affine" << std::endl;
    std::cout << "  --generate-config   Generate default config files and exit\n" << std::endl;
    std::cout << "  --sampling-percentage <0.0-1.0>  Sampling ratio (default 0.10)\n";
    
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << programName << " fixed.nrrd moving.nrrd output/" << std::endl;
    std::cout << "  " << programName << " --config Affine.json fixed.nrrd moving.nrrd output/" << std::endl;
    std::cout << "  " << programName << " --initial coarse.h5 fixed.nrrd moving.nrrd output/" << std::endl;
    std::cout << "  " << programName << " --transform Affine fixed.nrrd moving.nrrd output/\n" << std::endl;
    
    std::cout << "Supported image formats: NIFTI (.nii, .nii.gz), NRRD (.nrrd), MetaImage (.mhd/.mha)" << std::endl;
}

// ============================================================================
// 命令行参数解析
// ============================================================================

bool ParseCommandLine(int argc, std::vector<std::string>& args, CommandLineArgs& parsedArgs)
{
    std::vector<std::string> positionalArgs;
    
    for (size_t i = 1; i < args.size(); ++i)
    {
        const std::string& arg = args[i];
        
        if (arg == "--help" || arg == "-h")
        {
            parsedArgs.showHelp = true;
            return true;
        }
        else if (arg == "--generate-config")
        {
            parsedArgs.generateConfig = true;
            return true;
        }
        else if (arg == "--config")
        {
            if (i + 1 < args.size())
            {
                parsedArgs.configFilePath = args[++i];
            }
            else
            {
                std::cerr << "[Error] --config requires a file path" << std::endl;
                return false;
            }
        }
        else if (arg == "--initial")
        {
            if (i + 1 < args.size())
            {
                parsedArgs.initialTransformPath = args[++i];
            }
            else
            {
                std::cerr << "[Error] --initial requires a file path" << std::endl;
                return false;
            }
        }
        else if (arg == "--transform")
        {
            if (i + 1 < args.size())
            {
                parsedArgs.transformType = args[++i];
            }
            else
            {
                std::cerr << "[Error] --transform requires a type (Rigid or Affine)" << std::endl;
                return false;
            }
        }
        else if (arg == "--sampling-percentage")
        {
            if (i + 1 < args.size())
            {
                parsedArgs.samplingPercentage = std::stod(args[++i]);
            }
            else
            {
                std::cerr << "[Error] --sampling-percentage requires a numeric value (0.0-1.0)" << std::endl;
                return false;
            }
        }
        else if (arg == "--verbose")
        {
            parsedArgs.verbose = true;
        }
        else if (arg[0] == '-')
        {
            std::cerr << "[Error] Unknown option: " << arg << std::endl;
            return false;
        }
        else
        {
            positionalArgs.push_back(arg);
        }
    }
    
    if (positionalArgs.size() >= 3)
    {
        parsedArgs.fixedImagePath = positionalArgs[0];
        parsedArgs.movingImagePath = positionalArgs[1];
        parsedArgs.outputFolder = positionalArgs[2];
        return true;
    }
    else if (!parsedArgs.showHelp && !parsedArgs.generateConfig)
    {
        std::cerr << "[Error] Missing required arguments" << std::endl;
        return false;
    }
    
    return true;
}

// ============================================================================
// 主函数
// ============================================================================

int main(int argc, char* argv[])
{
    std::cout << "=== Mutual Information Registration (Custom Implementation) ===" << std::endl;
    std::cout << "Multi-Resolution 3D Registration for MRI/CT Images" << std::endl;
    std::cout << "============================================================\n" << std::endl;

    // 解析命令行参数
    std::vector<std::string> args;
    
#ifdef _WIN32
    int wargc;
    LPWSTR* wargv = CommandLineToArgvW(GetCommandLineW(), &wargc);
    for (int i = 0; i < wargc; ++i)
    {
        args.push_back(WideToUtf8(wargv[i]));
    }
    LocalFree(wargv);
#else
    for (int i = 0; i < argc; ++i)
    {
        args.push_back(argv[i]);
    }
#endif

    CommandLineArgs parsedArgs;
    if (!ParseCommandLine(argc, args, parsedArgs))
    {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    if (parsedArgs.showHelp)
    {
        PrintUsage(argv[0]);
        return EXIT_SUCCESS;
    }
    
    if (parsedArgs.generateConfig)
    {
        std::cout << "[Generating default configuration files...]" << std::endl;
        ConfigManager::CreateDefaultConfigFile("Rigid.json", ConfigManager::TransformType::Rigid);
        ConfigManager::CreateDefaultConfigFile("Affine.json", ConfigManager::TransformType::Affine);
        return EXIT_SUCCESS;
    }

    std::cout << "[Input Configuration]" << std::endl;
    std::cout << "  Fixed Image:  " << parsedArgs.fixedImagePath << std::endl;
    std::cout << "  Moving Image: " << parsedArgs.movingImagePath << std::endl;
    std::cout << "  Output Folder: " << parsedArgs.outputFolder << std::endl;
    if (!parsedArgs.configFilePath.empty())
    {
        std::cout << "  Config File: " << parsedArgs.configFilePath << std::endl;
    }
    if (!parsedArgs.initialTransformPath.empty())
    {
        std::cout << "  Initial Transform: " << parsedArgs.initialTransformPath << std::endl;
    }
    std::cout << "  Transform Type: " << parsedArgs.transformType << std::endl;
    std::cout << std::endl;

    // 验证输入文件
    if (!fs::exists(parsedArgs.fixedImagePath))
    {
        std::cerr << "[Error] Fixed image not found: " << parsedArgs.fixedImagePath << std::endl;
        return EXIT_FAILURE;
    }
    if (!fs::exists(parsedArgs.movingImagePath))
    {
        std::cerr << "[Error] Moving image not found: " << parsedArgs.movingImagePath << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        // 加载配置
        ConfigManager configManager;
        if (!parsedArgs.configFilePath.empty())
        {
            configManager.LoadFromFile(parsedArgs.configFilePath);
        }
        
        // 命令行覆盖配置文件
        if (!parsedArgs.transformType.empty())
        {
            configManager.SetTransformType(parsedArgs.transformType);
        }
        
        // 打印配置
        configManager.PrintConfig();
        
        // 创建配准对象
        ImageRegistration registration;
        
        // 从配置加载参数
        registration.LoadFromConfig(configManager.GetConfig());
        registration.SetVerbose(parsedArgs.verbose);

        // 命令行覆盖采样百分比
        if (parsedArgs.samplingPercentage >= 0.0)
        {
            registration.SetSamplingPercentage(parsedArgs.samplingPercentage);
        }
        
        // 加载图像
        std::cout << "\n[Loading Images...]" << std::endl;
        registration.SetFixedImagePath(parsedArgs.fixedImagePath);
        registration.SetMovingImagePath(parsedArgs.movingImagePath);
        
        // 加载初始变换(如果有)
        if (!parsedArgs.initialTransformPath.empty())
        {
            if (fs::exists(parsedArgs.initialTransformPath))
            {
                registration.LoadInitialTransform(parsedArgs.initialTransformPath);
            }
            else
            {
                std::cerr << "[Warning] Initial transform file not found: " 
                          << parsedArgs.initialTransformPath << std::endl;
            }
        }
        
        // 设置观察者
        registration.SetIterationObserver([](int iteration, double value, double stepLength) {
            std::cout << "  Iteration " << std::setw(4) << iteration 
                      << " | MI Value: " << std::fixed << std::setprecision(6) << value
                      << " | Step: " << std::scientific << std::setprecision(2) << stepLength
                      << std::endl;
        });
        
        registration.SetLevelObserver([](int level, int shrinkFactor, double sigma) {
            std::cout << "\n[Multi-Resolution Level " << (level + 1) << "]" << std::endl;
            std::cout << "  Shrink Factor: " << shrinkFactor << "x" << std::endl;
            std::cout << "  Smoothing Sigma: " << sigma << " mm" << std::endl;
        });
        
        // 执行配准
        std::cout << "\n[Starting Registration...]" << std::endl;
        registration.Update();
        
        // 获取结果
        double elapsedTime = registration.GetElapsedTime();
        
        std::cout << "\n[Registration Completed]" << std::endl;
        std::cout << "  Total Time: " << std::fixed << std::setprecision(2) 
                  << elapsedTime << " seconds" << std::endl;
        
        // 输出变换参数
        std::cout << "\n[Final Transform Parameters]" << std::endl;
        
        if (configManager.GetConfig().transformType == ConfigManager::TransformType::Rigid)
        {
            auto rigidTransform = registration.GetRigidTransform();
            auto parameters = rigidTransform->GetParameters();
            
            std::cout << "  Rotation (rad):    [" 
                      << std::fixed << std::setprecision(6)
                      << parameters[0] << ", " 
                      << parameters[1] << ", " 
                      << parameters[2] << "]" << std::endl;
            std::cout << "  Translation (mm):  [" 
                      << std::fixed << std::setprecision(4)
                      << parameters[3] << ", " 
                      << parameters[4] << ", " 
                      << parameters[5] << "]" << std::endl;
            
            auto center = rigidTransform->GetCenter();
            std::cout << "  Rotation Center:   [" 
                      << std::fixed << std::setprecision(2)
                      << center[0] << ", " 
                      << center[1] << ", " 
                      << center[2] << "]" << std::endl;
        }
        else
        {
            auto affineTransform = registration.GetAffineTransform();
            auto parameters = affineTransform->GetParameters();
            
            std::cout << "  Matrix:" << std::endl;
            std::cout << "    [" << std::fixed << std::setprecision(6)
                      << parameters[0] << ", " << parameters[1] << ", " << parameters[2] << "]" << std::endl;
            std::cout << "    [" << parameters[3] << ", " << parameters[4] << ", " << parameters[5] << "]" << std::endl;
            std::cout << "    [" << parameters[6] << ", " << parameters[7] << ", " << parameters[8] << "]" << std::endl;
            std::cout << "  Translation (mm):  [" 
                      << std::fixed << std::setprecision(4)
                      << parameters[9] << ", " 
                      << parameters[10] << ", " 
                      << parameters[11] << "]" << std::endl;
            
            auto center = affineTransform->GetCenter();
            std::cout << "  Center:   [" 
                      << std::fixed << std::setprecision(2)
                      << center[0] << ", " 
                      << center[1] << ", " 
                      << center[2] << "]" << std::endl;
        }
        
        // 保存变换
        std::cout << "\n[Saving Transform...]" << std::endl;
        
        // 关键修复：直接保存单个最终变换，不使用 CompositeTransform
        // 这避免了 .h5 文件中包含多个变换导致的混淆
        itk::Transform<double, 3, 3>::Pointer finalTransform;
        if (parsedArgs.transformType == "Rigid")
        {
            finalTransform = registration.GetRigidTransform();
        }
        else
        {
            finalTransform = registration.GetAffineTransform();
        }
        
        try
        {
            fs::path outputPath = fs::path(parsedArgs.outputFolder) / GenerateTimestampFilename();
            
            if (!fs::exists(parsedArgs.outputFolder))
            {
                fs::create_directories(parsedArgs.outputFolder);
            }
            
            using WriterType = itk::TransformFileWriter;
            auto writer = WriterType::New();
            writer->SetFileName(outputPath.string());
            writer->SetInput(finalTransform);
            writer->Update();
            
            std::cout << "[Transform Saved] " << outputPath.string() << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cerr << "[Error] Failed to save transform: " << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        
        std::cout << "\n=== Registration Successfully Completed ===" << std::endl;
        return EXIT_SUCCESS;
    }
    catch (const itk::ExceptionObject& e)
    {
        std::cerr << "\n[ITK Error] " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::exception& e)
    {
        std::cerr << "\n[Error] " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}

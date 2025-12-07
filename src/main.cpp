#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

#include "ImageRegistration.h"

#include <itkTransformFileWriter.h>

namespace fs = std::filesystem;

// 函数声明
std::string GenerateTimestampFilename();
bool SaveTransformAsH5(ImageRegistration::TransformPointer transform, const std::string& outputFolder);
void PrintUsage(const char* programName);

// Windows 宽字符转 UTF-8
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

// 生成带时间戳的输出文件名
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

// 保存变换到 HDF5 格式
bool SaveTransformAsH5(ImageRegistration::TransformPointer transform, 
                       const std::string& outputFolder)
{
    try
    {
        fs::path outputPath = fs::path(outputFolder) / GenerateTimestampFilename();
        
        // 确保输出文件夹存在
        if (!fs::exists(outputFolder))
        {
            fs::create_directories(outputFolder);
        }
        
        // 使用 ITK 保存变换
        using WriterType = itk::TransformFileWriter;
        auto writer = WriterType::New();
        writer->SetFileName(outputPath.string());
        writer->SetInput(transform);
        writer->Update();
        
        std::cout << "[Transform Saved] " << outputPath.string() << std::endl;
        return true;
    }
    catch (const std::exception& e)
    {
        std::cerr << "[Error] Failed to save transform: " << e.what() << std::endl;
        return false;
    }
}

void PrintUsage(const char* programName)
{
    std::cout << "Usage: " << programName << " <fixed_image> <moving_image> <output_folder>\n"
              << "\n"
              << "Arguments:\n"
              << "  fixed_image   - Path to the fixed (reference) image (MRI/CT)\n"
              << "  moving_image  - Path to the moving image to be registered\n"
              << "  output_folder - Folder path for output .h5 transform file\n"
              << "\n"
              << "Supported formats: NIFTI (.nii, .nii.gz), NRRD (.nrrd), MetaImage (.mhd/.mha)\n"
              << std::endl;
}

int main(int argc, char* argv[])
{
    std::cout << "=== Mutual Information Registration (Custom Implementation) ===" << std::endl;
    std::cout << "Multi-Resolution 3D Rigid Registration for MRI/CT Images" << std::endl;
    std::cout << "============================================================\n" << std::endl;

#ifdef _WIN32
    // Windows: 使用宽字符 API 支持中文路径
    int wargc;
    LPWSTR* wargv = CommandLineToArgvW(GetCommandLineW(), &wargc);
    
    if (wargc < 4)
    {
        LocalFree(wargv);
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    std::string fixedImagePath = WideToUtf8(wargv[1]);
    std::string movingImagePath = WideToUtf8(wargv[2]);
    std::string outputFolder = WideToUtf8(wargv[3]);
    LocalFree(wargv);
#else
    if (argc < 4)
    {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    std::string fixedImagePath = argv[1];
    std::string movingImagePath = argv[2];
    std::string outputFolder = argv[3];
#endif

    std::cout << "[Input Configuration]" << std::endl;
    std::cout << "  Fixed Image:  " << fixedImagePath << std::endl;
    std::cout << "  Moving Image: " << movingImagePath << std::endl;
    std::cout << "  Output Folder: " << outputFolder << std::endl;
    std::cout << std::endl;

    // 验证输入文件存在
    if (!fs::exists(fixedImagePath))
    {
        std::cerr << "[Error] Fixed image not found: " << fixedImagePath << std::endl;
        return EXIT_FAILURE;
    }
    if (!fs::exists(movingImagePath))
    {
        std::cerr << "[Error] Moving image not found: " << movingImagePath << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        // 创建配准对象
        ImageRegistration registration;
        
        // 设置图像
        std::cout << "[Loading Images...]" << std::endl;
        registration.SetFixedImagePath(fixedImagePath);
        registration.SetMovingImagePath(movingImagePath);
        
        // 配置多分辨率参数
        registration.SetNumberOfLevels(3);
        registration.SetShrinkFactors({4, 2, 1});
        registration.SetSmoothingSigmas({2.0, 1.0, 0.0});
        
        // 配置优化器参数
        registration.SetLearningRate(0.5);
        registration.SetMinimumStepLength(0.0001);
        registration.SetNumberOfIterations(300);
        registration.SetRelaxationFactor(0.8);
        registration.SetGradientMagnitudeTolerance(1e-4);
        
        // 配置度量参数
        registration.SetNumberOfHistogramBins(50);
        registration.SetNumberOfSamples(100000);
        
        // 设置迭代观察者
        registration.SetIterationObserver([](int iteration, double value, double stepLength) {
            std::cout << "  Iteration " << std::setw(4) << iteration 
                      << " | MI Value: " << std::fixed << std::setprecision(6) << value
                      << " | Step: " << std::scientific << std::setprecision(2) << stepLength
                      << std::endl;
        });
        
        // 设置多分辨率级别观察者
        registration.SetLevelObserver([](int level, int shrinkFactor, double sigma) {
            std::cout << "\n[Multi-Resolution Level " << (level + 1) << "]" << std::endl;
            std::cout << "  Shrink Factor: " << shrinkFactor << "x" << std::endl;
            std::cout << "  Smoothing Sigma: " << sigma << " mm" << std::endl;
        });
        
        // 执行配准
        std::cout << "\n[Starting Registration...]" << std::endl;
        registration.Update();
        
        // 获取结果
        auto finalTransform = registration.GetFinalTransform();
        double elapsedTime = registration.GetElapsedTime();
        
        // 输出结果
        std::cout << "\n[Registration Completed]" << std::endl;
        std::cout << "  Total Time: " << std::fixed << std::setprecision(2) 
                  << elapsedTime << " seconds" << std::endl;
        
        // 输出变换参数
        auto parameters = finalTransform->GetParameters();
        std::cout << "\n[Final Transform Parameters]" << std::endl;
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
        
        auto center = finalTransform->GetCenter();
        std::cout << "  Rotation Center:   [" 
                  << std::fixed << std::setprecision(2)
                  << center[0] << ", " 
                  << center[1] << ", " 
                  << center[2] << "]" << std::endl;
        
        // 保存变换
        std::cout << "\n[Saving Transform...]" << std::endl;
        if (!SaveTransformAsH5(finalTransform, outputFolder))
        {
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

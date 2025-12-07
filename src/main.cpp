#include <iostream>
#include <string>
#include <sstream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <ctime>

#ifdef _WIN32
#include <windows.h>
#endif

#include <itkImageFileReader.h>
#include <itkImage.h>
#include <itkEuler3DTransform.h>
#include <itkRegularStepGradientDescentOptimizerv4.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
#include <itkImageRegistrationMethodv4.h>
#include <itkCenteredTransformInitializer.h>
#include <itkTransformFileWriter.h>

using ImageType = itk::Image<float, 3>;
using TransformType = itk::Euler3DTransform<double>;

// Forward declaration
bool RunMattesRigidRegistration(const std::string &fixedFile, const std::string &movingFile, 
                                TransformType::Pointer &outTransform, double &elapsedSeconds);

#ifdef _WIN32
static std::string WideToUtf8(const std::wstring &wstr)
{
    if (wstr.empty()) return std::string();
    int size_needed = WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), -1, NULL, 0, NULL, NULL);
    std::string strTo(size_needed, '\0');
    WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), -1, &strTo[0], size_needed, NULL, NULL);
    if (!strTo.empty() && strTo.back() == '\0') strTo.pop_back();
    return strTo;
}
#endif

// Save transform as HDF5 format
bool SaveTransformAsH5(TransformType::Pointer transform, const std::string &h5FilePath)
{
    try {
        auto writer = itk::TransformFileWriterTemplate<double>::New();
        writer->SetFileName(h5FilePath);
        writer->SetInput(transform);
        writer->Update();
        return true;
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "Error saving transform to H5: " << ex << std::endl;
        return false;
    }
}

#ifdef _WIN32
int mainImplW(int argc, wchar_t** argvW)
{
    namespace fs = std::filesystem;

    if (argc < 4) {
        std::cout << "Usage: MIRegistration <fixed_image> <moving_image> <output_folder>\n";
        std::cout << "  Output will be saved as transform_YYYYMMDD_HHMMSS.h5 in output_folder\n";
        return 0;
    }

    fs::path fixedPath(argvW[1]);
    fs::path movingPath(argvW[2]);
    fs::path outputFolder(argvW[3]);

    std::cout << "\n=== Mutual Information Registration ===\n";
    std::cout << "Fixed:  " << WideToUtf8(fixedPath.wstring()) << std::endl;
    std::cout << "Moving: " << WideToUtf8(movingPath.wstring()) << std::endl;
    std::cout << "Output: " << WideToUtf8(outputFolder.wstring()) << std::endl;

    if (!fs::exists(fixedPath)) {
        std::cerr << "ERROR: Fixed image file not found: " << WideToUtf8(fixedPath.wstring()) << std::endl;
        return 1;
    }
    if (!fs::exists(movingPath)) {
        std::cerr << "ERROR: Moving image file not found: " << WideToUtf8(movingPath.wstring()) << std::endl;
        return 1;
    }

    // Create output folder if not exists
    if (!fs::exists(outputFolder)) {
        fs::create_directories(outputFolder);
    }

    // Verify images can be loaded and print metadata
    ImageType::Pointer fixedImage, movingImage;
    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(fixedPath.u8string());
        reader->Update();
        fixedImage = reader->GetOutput();
        std::cout << "Fixed image loaded successfully.\n";
        
        auto origin = fixedImage->GetOrigin();
        auto spacing = fixedImage->GetSpacing();
        auto size = fixedImage->GetLargestPossibleRegion().GetSize();
        auto direction = fixedImage->GetDirection();
        
        std::cout << "  Size: " << size[0] << " x " << size[1] << " x " << size[2] << "\n";
        std::cout << "  Spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << " mm\n";
        std::cout << "  Origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << " mm\n";
        std::cout << "  Direction Matrix:\n";
        for (unsigned i = 0; i < 3; ++i) {
            std::cout << "    " << direction[i][0] << ", " << direction[i][1] << ", " << direction[i][2] << "\n";
        }
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "ERROR loading fixed image:\n" << ex << std::endl;
        return 1;
    }

    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(movingPath.u8string());
        reader->Update();
        movingImage = reader->GetOutput();
        std::cout << "\nMoving image loaded successfully.\n";
        
        auto origin = movingImage->GetOrigin();
        auto spacing = movingImage->GetSpacing();
        auto size = movingImage->GetLargestPossibleRegion().GetSize();
        auto direction = movingImage->GetDirection();
        
        std::cout << "  Size: " << size[0] << " x " << size[1] << " x " << size[2] << "\n";
        std::cout << "  Spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << " mm\n";
        std::cout << "  Origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << " mm\n";
        std::cout << "  Direction Matrix:\n";
        for (unsigned i = 0; i < 3; ++i) {
            std::cout << "    " << direction[i][0] << ", " << direction[i][1] << ", " << direction[i][2] << "\n";
        }
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "ERROR loading moving image:\n" << ex << std::endl;
        return 1;
    }

    // Run registration
    std::cout << "\n--- Starting Registration ---\n";
    TransformType::Pointer transform;
    double elapsed = 0.0;
    if (!RunMattesRigidRegistration(WideToUtf8(fixedPath.wstring()), 
                                     WideToUtf8(movingPath.wstring()), 
                                     transform, elapsed)) {
        std::cerr << "ERROR: Registration failed.\n";
        return 1;
    }

    // Print transform parameters and 4x4 matrix
    std::cout << "\n=== Registration Results ===\n";
    std::cout << "Time elapsed: " << std::fixed << std::setprecision(2) << elapsed << " seconds\n";
    
    // Print detailed transform parameters
    std::cout << "\nTransform Parameters:\n";
    std::cout << "  Rotation Center: " 
              << transform->GetCenter()[0] << ", " 
              << transform->GetCenter()[1] << ", " 
              << transform->GetCenter()[2] << " mm\n";
    std::cout << "  Rotation (rad): " 
              << transform->GetAngleX() << ", " 
              << transform->GetAngleY() << ", " 
              << transform->GetAngleZ() << "\n";
    std::cout << "  Translation: " 
              << transform->GetTranslation()[0] << ", " 
              << transform->GetTranslation()[1] << ", " 
              << transform->GetTranslation()[2] << " mm\n";
    
    // Compute and print full 4x4 transformation matrix
    // This matrix includes the rotation center offset
    auto matrix = transform->GetMatrix();
    auto offset = transform->GetOffset();
    std::cout << "\n4x4 Transformation Matrix (with rotation center):\n";
    for (unsigned r = 0; r < 3; ++r) {
        std::cout << std::setw(12) << std::setprecision(6) << matrix[r][0] << " "
                  << std::setw(12) << matrix[r][1] << " "
                  << std::setw(12) << matrix[r][2] << " "
                  << std::setw(12) << offset[r] << "\n";
    }
    std::cout << std::setw(12) << 0.0 << " "
              << std::setw(12) << 0.0 << " "
              << std::setw(12) << 0.0 << " "
              << std::setw(12) << 1.0 << "\n";

    // Generate output filename with timestamp
    auto now = std::chrono::system_clock::now();
    auto now_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
    localtime_s(&tm_buf, &now_t);
    
    std::ostringstream oss;
    oss << "transform_" 
        << std::setfill('0') << std::setw(4) << (tm_buf.tm_year + 1900)
        << std::setw(2) << (tm_buf.tm_mon + 1)
        << std::setw(2) << tm_buf.tm_mday << "_"
        << std::setw(2) << tm_buf.tm_hour
        << std::setw(2) << tm_buf.tm_min
        << std::setw(2) << tm_buf.tm_sec
        << ".h5";
    
    fs::path h5Path = outputFolder / oss.str();
    
    // Save transform as HDF5
    std::cout << "\nSaving transform to: " << WideToUtf8(h5Path.wstring()) << std::endl;
    if (SaveTransformAsH5(transform, h5Path.u8string())) {
        std::cout << "Transform saved successfully!\n";
    } else {
        std::cerr << "ERROR: Failed to save transform.\n";
        return 1;
    }

    return 0;
}
#endif

int mainImpl(int argc, char** argv) {
    namespace fs = std::filesystem;

    if (argc < 4) {
        std::cout << "Usage: MIRegistration <fixed_image> <moving_image> <output_folder>\n";
        std::cout << "  Output will be saved as transform_YYYYMMDD_HHMMSS.h5 in output_folder\n";
        return 0;
    }

    fs::path fixedPath = argv[1];
    fs::path movingPath = argv[2];
    fs::path outputFolder = argv[3];

    std::cout << "\n=== Mutual Information Registration ===\n";
    std::cout << "Fixed:  " << fixedPath.string() << std::endl;
    std::cout << "Moving: " << movingPath.string() << std::endl;
    std::cout << "Output: " << outputFolder.string() << std::endl;

    if (!fs::exists(fixedPath)) {
        std::cerr << "ERROR: Fixed image file not found: " << fixedPath.string() << std::endl;
        return 1;
    }
    if (!fs::exists(movingPath)) {
        std::cerr << "ERROR: Moving image file not found: " << movingPath.string() << std::endl;
        return 1;
    }

    // Create output folder if not exists
    if (!fs::exists(outputFolder)) {
        fs::create_directories(outputFolder);
    }

    // Verify images can be loaded and print metadata
    ImageType::Pointer fixedImage, movingImage;
    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(fixedPath.string());
        reader->Update();
        fixedImage = reader->GetOutput();
        std::cout << "Fixed image loaded successfully.\n";
        
        auto origin = fixedImage->GetOrigin();
        auto spacing = fixedImage->GetSpacing();
        auto size = fixedImage->GetLargestPossibleRegion().GetSize();
        auto direction = fixedImage->GetDirection();
        
        std::cout << "  Size: " << size[0] << " x " << size[1] << " x " << size[2] << "\n";
        std::cout << "  Spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << " mm\n";
        std::cout << "  Origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << " mm\n";
        std::cout << "  Direction Matrix:\n";
        for (unsigned i = 0; i < 3; ++i) {
            std::cout << "    " << direction[i][0] << ", " << direction[i][1] << ", " << direction[i][2] << "\n";
        }
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "ERROR loading fixed image:\n" << ex << std::endl;
        return 1;
    }

    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(movingPath.string());
        reader->Update();
        movingImage = reader->GetOutput();
        std::cout << "\nMoving image loaded successfully.\n";
        
        auto origin = movingImage->GetOrigin();
        auto spacing = movingImage->GetSpacing();
        auto size = movingImage->GetLargestPossibleRegion().GetSize();
        auto direction = movingImage->GetDirection();
        
        std::cout << "  Size: " << size[0] << " x " << size[1] << " x " << size[2] << "\n";
        std::cout << "  Spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << " mm\n";
        std::cout << "  Origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << " mm\n";
        std::cout << "  Direction Matrix:\n";
        for (unsigned i = 0; i < 3; ++i) {
            std::cout << "    " << direction[i][0] << ", " << direction[i][1] << ", " << direction[i][2] << "\n";
        }
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "ERROR loading moving image:\n" << ex << std::endl;
        return 1;
    }

    // Run registration
    std::cout << "\n--- Starting Registration ---\n";
    TransformType::Pointer transform;
    double elapsed = 0.0;
    if (!RunMattesRigidRegistration(fixedPath.string(), movingPath.string(), transform, elapsed)) {
        std::cerr << "ERROR: Registration failed.\n";
        return 1;
    }

    // Print transform parameters and 4x4 matrix
    std::cout << "\n=== Registration Results ===\n";
    std::cout << "Time elapsed: " << std::fixed << std::setprecision(2) << elapsed << " seconds\n";
    
    // Print detailed transform parameters
    std::cout << "\nTransform Parameters:\n";
    std::cout << "  Rotation Center: " 
              << transform->GetCenter()[0] << ", " 
              << transform->GetCenter()[1] << ", " 
              << transform->GetCenter()[2] << " mm\n";
    std::cout << "  Rotation (rad): " 
              << transform->GetAngleX() << ", " 
              << transform->GetAngleY() << ", " 
              << transform->GetAngleZ() << "\n";
    std::cout << "  Translation: " 
              << transform->GetTranslation()[0] << ", " 
              << transform->GetTranslation()[1] << ", " 
              << transform->GetTranslation()[2] << " mm\n";
    
    // Compute and print full 4x4 transformation matrix
    // This matrix includes the rotation center offset
    auto matrix = transform->GetMatrix();
    auto offset = transform->GetOffset();
    std::cout << "\n4x4 Transformation Matrix (with rotation center):\n";
    for (unsigned r = 0; r < 3; ++r) {
        std::cout << std::setw(12) << std::setprecision(6) << matrix[r][0] << " "
                  << std::setw(12) << matrix[r][1] << " "
                  << std::setw(12) << matrix[r][2] << " "
                  << std::setw(12) << offset[r] << "\n";
    }
    std::cout << std::setw(12) << 0.0 << " "
              << std::setw(12) << 0.0 << " "
              << std::setw(12) << 0.0 << " "
              << std::setw(12) << 1.0 << "\n";

    // Generate output filename with timestamp
    auto now = std::chrono::system_clock::now();
    auto now_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
#ifdef _WIN32
    localtime_s(&tm_buf, &now_t);
#else
    localtime_r(&now_t, &tm_buf);
#endif
    
    std::ostringstream oss;
    oss << "transform_" 
        << std::setfill('0') << std::setw(4) << (tm_buf.tm_year + 1900)
        << std::setw(2) << (tm_buf.tm_mon + 1)
        << std::setw(2) << tm_buf.tm_mday << "_"
        << std::setw(2) << tm_buf.tm_hour
        << std::setw(2) << tm_buf.tm_min
        << std::setw(2) << tm_buf.tm_sec
        << ".h5";
    
    fs::path h5Path = outputFolder / oss.str();
    
    // Save transform as HDF5
    std::cout << "\nSaving transform to: " << h5Path.string() << std::endl;
    if (SaveTransformAsH5(transform, h5Path.string())) {
        std::cout << "Transform saved successfully!\n";
    } else {
        std::cerr << "ERROR: Failed to save transform.\n";
        return 1;
    }

    return 0;
}

int main(int argc, char** argv) {
#ifdef _WIN32
    int wargc = 0;
    wchar_t **wargv = CommandLineToArgvW(GetCommandLineW(), &wargc);
    if (wargv != nullptr) {
        int ret = mainImplW(wargc, wargv);
        LocalFree(wargv);
        return ret;
    }
#endif
    return mainImpl(argc, argv);
}

bool RunMattesRigidRegistration(const std::string &fixedFile, const std::string &movingFile, TransformType::Pointer &outTransform, double &elapsedSeconds)
{
    using FixedImageType = ImageType;
    using MovingImageType = ImageType;

    // Readers
    auto fixedReader = itk::ImageFileReader<FixedImageType>::New();
    fixedReader->SetFileName(fixedFile);
    auto movingReader = itk::ImageFileReader<MovingImageType>::New();
    movingReader->SetFileName(movingFile);

    try {
        fixedReader->Update();
        movingReader->Update();
    }
    catch (const itk::ExceptionObject &ex) {
        std::cerr << "Image read error: " << ex << std::endl;
        return false;
    }

    auto fixedImage = fixedReader->GetOutput();
    auto movingImage = movingReader->GetOutput();

    // Initial transform (centered)
    auto initialTransform = TransformType::New();
    using InitializerType = itk::CenteredTransformInitializer<TransformType, FixedImageType, MovingImageType>;
    auto initializer = InitializerType::New();
    initializer->SetTransform(initialTransform);
    initializer->SetFixedImage(fixedImage);
    initializer->SetMovingImage(movingImage);
    initializer->GeometryOn();  // 使用几何中心而不是质心,更稳定
    initializer->InitializeTransform();

    std::cout << "Initial transform center: " << initialTransform->GetCenter() << std::endl;

    // Metric
    using MetricType = itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType>;
    auto metric = MetricType::New();
    metric->SetNumberOfHistogramBins(50);
    metric->SetUseMovingImageGradientFilter(false);  // 提高稳定性
    metric->SetUseFixedImageGradientFilter(false);

    // Optimizer
    using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
    auto optimizer = OptimizerType::New();
    optimizer->SetLearningRate(0.5);  // 降低学习率提高稳定性
    optimizer->SetMinimumStepLength(0.0001);
    optimizer->SetNumberOfIterations(300);  // 增加迭代次数
    optimizer->SetRelaxationFactor(0.8);  // 增加松弛因子
    optimizer->SetGradientMagnitudeTolerance(1e-4);
    optimizer->SetReturnBestParametersAndValue(true);  // 返回最佳参数

    // Registration - 需要先定义,因为Observer会引用它
    using RegistrationType = itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, TransformType>;
    auto registration = RegistrationType::New();
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetInitialTransform(initialTransform);

    // Add observer to print iteration progress
    class IterationObserver : public itk::Command
    {
    public:
        using Self = IterationObserver;
        using Superclass = itk::Command;
        using Pointer = itk::SmartPointer<Self>;
        itkNewMacro(Self);

    protected:
        IterationObserver() = default;

    public:
        void Execute(itk::Object *caller, const itk::EventObject &event) override
        {
            Execute((const itk::Object *)caller, event);
        }

        void Execute(const itk::Object *object, const itk::EventObject &event) override
        {
            auto optimizer = static_cast<const OptimizerType *>(object);
            if (!itk::IterationEvent().CheckEvent(&event))
            {
                return;
            }

            unsigned int iter = optimizer->GetCurrentIteration();
            double value = optimizer->GetValue();
            double learningRate = optimizer->GetLearningRate();

            // Print progress every 10 iterations or at the first iteration
            if (iter % 10 == 0 || iter == 0)
            {
                std::cout << "  Iter: " << std::setw(4) << iter 
                          << "  Metric: " << std::setw(12) << std::fixed << std::setprecision(6) << value
                          << "  LearningRate: " << std::setw(10) << std::scientific << std::setprecision(4) << learningRate
                          << std::endl;
            }
        }
    };

    // Add observer for multi-resolution level changes
    class MultiResolutionObserver : public itk::Command
    {
    public:
        using Self = MultiResolutionObserver;
        using Superclass = itk::Command;
        using Pointer = itk::SmartPointer<Self>;
        itkNewMacro(Self);

    protected:
        MultiResolutionObserver() = default;

    public:
        void Execute(itk::Object *caller, const itk::EventObject &event) override
        {
            Execute((const itk::Object *)caller, event);
        }

        void Execute(const itk::Object *object, const itk::EventObject &event) override
        {
            if (!itk::MultiResolutionIterationEvent().CheckEvent(&event))
            {
                return;
            }

            auto registration = static_cast<const RegistrationType *>(object);
            unsigned int level = registration->GetCurrentLevel();
            unsigned int totalLevels = registration->GetNumberOfLevels();
            
            std::cout << "\n=== Multi-Resolution Level " << level << " of " << (totalLevels - 1) << " ===" << std::endl;
        }
    };

    auto iterObserver = IterationObserver::New();
    optimizer->AddObserver(itk::IterationEvent(), iterObserver);

    auto multiResObserver = MultiResolutionObserver::New();
    registration->AddObserver(itk::MultiResolutionIterationEvent(), multiResObserver);
    
    // Configure multi-resolution pyramid (默认是3层)
    // 你可以修改这些参数来调整金字塔策略
    constexpr unsigned int numberOfLevels = 3;
    RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize(numberOfLevels);
    shrinkFactorsPerLevel[0] = 4;  // 第一层: 下采样4倍 (粗略对齐)
    shrinkFactorsPerLevel[1] = 2;  // 第二层: 下采样2倍 (中等精度)
    shrinkFactorsPerLevel[2] = 1;  // 第三层: 原始分辨率 (精细对齐)
    
    RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize(numberOfLevels);
    smoothingSigmasPerLevel[0] = 2.0;  // 第一层: 高斯平滑2mm
    smoothingSigmasPerLevel[1] = 1.0;  // 第二层: 高斯平滑1mm
    smoothingSigmasPerLevel[2] = 0.0;  // 第三层: 不平滑
    
    registration->SetNumberOfLevels(numberOfLevels);
    registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
    registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
    registration->SmoothingSigmasAreSpecifiedInPhysicalUnitsOn();  // 平滑参数单位是mm

    std::cout << "\nMulti-Resolution Strategy: " << numberOfLevels << " levels" << std::endl;
    std::cout << "  Level 0: Shrink " << shrinkFactorsPerLevel[0] << "x, Smooth " << smoothingSigmasPerLevel[0] << " mm (coarse)" << std::endl;
    std::cout << "  Level 1: Shrink " << shrinkFactorsPerLevel[1] << "x, Smooth " << smoothingSigmasPerLevel[1] << " mm (medium)" << std::endl;
    std::cout << "  Level 2: Shrink " << shrinkFactorsPerLevel[2] << "x, Smooth " << smoothingSigmasPerLevel[2] << " mm (fine)" << std::endl;

    // 使用固定种子的确定性采样,确保可重复性
    auto fixedRegion = fixedImage->GetLargestPossibleRegion();
    auto fixedSize = fixedRegion.GetSize();
    unsigned long totalVoxels = fixedSize[0] * fixedSize[1] * fixedSize[2];
    unsigned long maxSamples = 100000;  // 增加采样数提高稳定性
    unsigned long numSamples = (totalVoxels < maxSamples) ? totalVoxels : maxSamples;
    double samplingPercentage = (double)numSamples / (double)totalVoxels;
    
    registration->SetMetricSamplingPercentage(samplingPercentage);
    registration->SetMetricSamplingStrategy(RegistrationType::RANDOM);
    registration->MetricSamplingReinitializeSeed(121212);  // 固定随机种子!

    std::cout << "Using " << numSamples << " samples (" 
              << std::fixed << std::setprecision(1) << (samplingPercentage * 100) 
              << "% of total voxels)\n";

    auto startTime = std::chrono::high_resolution_clock::now();

    try {
        registration->Update();
    }
    catch (const itk::ExceptionObject &ex) {
        std::cerr << "Registration error: " << ex << std::endl;
        return false;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    elapsedSeconds = std::chrono::duration<double>(endTime - startTime).count();

    // 获取最优值
    auto bestValue = optimizer->GetValue();
    std::cout << "Final metric value: " << bestValue << std::endl;

    outTransform = const_cast<TransformType*>(registration->GetOutput()->Get());

    return true;
}

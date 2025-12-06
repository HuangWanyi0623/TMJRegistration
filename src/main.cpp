#include "ImageRegistration.h"
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

#include <itkImageFileReader.h>
#include <itkImage.h>
#include <itkEuler3DTransform.h>
#include <itkRegularStepGradientDescentOptimizerv4.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
#include <itkImageRegistrationMethodv4.h>
#include <itkCenteredTransformInitializer.h>

using ImageType = itk::Image<float, 3>;
using TransformType = itk::Euler3DTransform<double>;

// Forward declaration
bool RunMattesRigidRegistration(const std::string &fixedFile, const std::string &movingFile, TransformType::Pointer &outTransform, double &elapsedSeconds);

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

int mainImplW(int argc, wchar_t** argvW)
{
    namespace fs = std::filesystem;

    if (argc < 4) {
        std::cout << "Usage: MIRegistration <fixed_image> <moving_image> <output_image>\n";
        return 0;
    }

    fs::path fixedPath(argvW[1]);
    fs::path movingPath(argvW[2]);
    fs::path outputPath(argvW[3]);

    std::cout << "Fixed:  " << WideToUtf8(fixedPath.wstring()) << std::endl;
    std::cout << "Moving: " << WideToUtf8(movingPath.wstring()) << std::endl;
    std::cout << "Output: " << WideToUtf8(outputPath.wstring()) << std::endl;

    if (!fs::exists(fixedPath)) {
        std::cerr << "Fixed image file not found: " << WideToUtf8(fixedPath.wstring()) << std::endl;
        return 1;
    }
    if (!fs::exists(movingPath)) {
        std::cerr << "Moving image file not found: " << WideToUtf8(movingPath.wstring()) << std::endl;
        return 1;
    }

    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(fixedPath.u8string());
        reader->Update();
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "Error loading fixed image: " << ex << std::endl;
        return 1;
    }

    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(movingPath.u8string());
        reader->Update();
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "Error loading moving image: " << ex << std::endl;
        return 1;
    }

    std::cout << "Images loaded successfully." << std::endl;

    TransformType::Pointer transform;
    double elapsed = 0.0;
    if (!RunMattesRigidRegistration(WideToUtf8(fixedPath.wstring()), WideToUtf8(movingPath.wstring()), transform, elapsed)) {
        std::cerr << "Registration failed." << std::endl;
        return 1;
    }

    // Print 4x4 matrix
    auto matrix = transform->GetMatrix();
    auto offset = transform->GetOffset();
    std::cout << "Registration completed successfully. Time: " << elapsed << " s" << std::endl;
    std::cout << "Transform (4x4):\n";
    for (unsigned r = 0; r < 3; ++r) {
        std::cout << matrix[r][0] << " " << matrix[r][1] << " " << matrix[r][2] << " " << offset[r] << std::endl;
    }
    std::cout << "0 0 0 1" << std::endl;

    return 0;
}
#endif

int mainImpl(int argc, char** argv) {
    if (argc < 4) {
        std::cout << "Usage: MIRegistration <fixed_image> <moving_image> <output_image>\n";
        return 0;
    }

    namespace fs = std::filesystem;

    fs::path fixed = argv[1];
    fs::path moving = argv[2];
    fs::path output = argv[3];

    std::cout << "Fixed:  " << fixed.string() << std::endl;
    std::cout << "Moving: " << moving.string() << std::endl;
    std::cout << "Output: " << output.string() << std::endl;

    if (!fs::exists(fixed)) {
        std::cerr << "Fixed image file not found: " << fixed.string() << std::endl;
        return 1;
    }
    if (!fs::exists(moving)) {
        std::cerr << "Moving image file not found: " << moving.string() << std::endl;
        return 1;
    }

    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(fixed.string());
        reader->Update();
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "Error loading fixed image: " << ex << std::endl;
        return 1;
    }

    try {
        auto reader = itk::ImageFileReader<ImageType>::New();
        reader->SetFileName(moving.string());
        reader->Update();
    }
    catch (const itk::ExceptionObject& ex) {
        std::cerr << "Error loading moving image: " << ex << std::endl;
        return 1;
    }

    std::cout << "Images loaded successfully." << std::endl;

    TransformType::Pointer transform;
    double elapsed = 0.0;
    if (!RunMattesRigidRegistration(fixed.string(), moving.string(), transform, elapsed))
    {
        std::cerr << "Registration failed." << std::endl;
        return 1;
    }

    auto matrix = transform->GetMatrix();
    auto offset = transform->GetOffset();
    std::cout << "Registration completed successfully. Time: " << elapsed << " s" << std::endl;
    std::cout << "Transform (4x4):\n";
    for (unsigned r = 0; r < 3; ++r) {
        std::cout << matrix[r][0] << " " << matrix[r][1] << " " << matrix[r][2] << " " << offset[r] << std::endl;
    }
    std::cout << "0 0 0 1" << std::endl;

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
    auto transform = TransformType::New();
    transform->SetIdentity();

    using TransformInitializerType = itk::CenteredTransformInitializer<TransformType, FixedImageType, MovingImageType>;
    auto initializer = TransformInitializerType::New();
    initializer->SetTransform(transform);
    initializer->SetFixedImage(fixedImage);
    initializer->SetMovingImage(movingImage);
    initializer->MomentsOn();
    initializer->InitializeTransform();

    // Metric
    using MetricType = itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType>;
    auto metric = MetricType::New();
    metric->SetNumberOfHistogramBins(32);

    // Set number of samples adaptively based on image size (use percentage sampling on registration)
    auto regionSize = fixedImage->GetLargestPossibleRegion().GetSize();
    unsigned long long totalVoxels = static_cast<unsigned long long>(regionSize[0]) * regionSize[1] * regionSize[2];
    unsigned long long spatialSamples = std::min<unsigned long long>(totalVoxels, 200000ULL);
    double samplingPercentage = static_cast<double>(spatialSamples) / static_cast<double>(totalVoxels);

    // Optimizer
    using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
    auto optimizer = OptimizerType::New();
    // Increase maximum iterations for a more substantial optimization
    optimizer->SetNumberOfIterations(500);
    optimizer->SetMinimumStepLength(1e-4);
    // removed SetMaximumStepLength (not available in this optimizer class)

    // Registration
    using RegistrationType = itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, TransformType>;
    auto registration = RegistrationType::New();
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);
    registration->SetInitialTransform(transform);
    registration->InPlaceOn();

    // Register samplings
    registration->SetMetricSamplingStrategy(RegistrationType::RANDOM);
    registration->SetMetricSamplingPercentage(samplingPercentage);

    try {
        auto tStart = std::chrono::steady_clock::now();
        registration->Update();
        auto tEnd = std::chrono::steady_clock::now();
        elapsedSeconds = std::chrono::duration<double>(tEnd - tStart).count();
    }
    catch (const itk::ExceptionObject &ex) {
        std::cerr << "Registration failed: " << ex << std::endl;
        return false;
    }

    // Get final parameters and set transform
    auto finalParams = optimizer->GetCurrentPosition();
    transform->SetParameters(finalParams);

    outTransform = transform;

    return true;
}

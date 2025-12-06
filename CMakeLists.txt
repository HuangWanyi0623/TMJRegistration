cmake_minimum_required(VERSION 3.16)
project(MIRegistration)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 使用vcpkg工具链
# 运行cmake时需要指定: -DCMAKE_TOOLCHAIN_FILE=[vcpkg路径]/scripts/buildsystems/vcpkg.cmake

# 查找ITK (通过vcpkg安装)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# 包含头文件目录
include_directories(${CMAKE_SOURCE_DIR}/include)

# 源文件
set(SOURCES
    src/MattesMutualInformation.cpp
    src/RegularStepGradientDescentOptimizer.cpp
    src/ImageRegistration.cpp
    src/main.cpp
)

# 头文件
set(HEADERS
    include/MattesMutualInformation.h
    include/RegularStepGradientDescentOptimizer.h
    include/ImageRegistration.h
)

# 创建可执行文件
add_executable(MIRegistration ${SOURCES} ${HEADERS})

# 针对 MSVC 明确指定源文件编码为 UTF-8，避免C4819和字符串相关错误
if(MSVC)
    target_compile_options(MIRegistration PRIVATE "/utf-8")
endif()

# 链接ITK库
target_link_libraries(MIRegistration ${ITK_LIBRARIES})

# 设置输出目录
set_target_properties(MIRegistration PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

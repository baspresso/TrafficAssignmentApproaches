cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER "C:/mingw64/bin/gcc.exe")

project(TrafficAssignmentApproaches)

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_HOST_WIN32)
    set(CMAKE_TOOLCHAIN_FILE "C:/tools/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
endif()

set(EIGEN3_INCLUDE_DIR "C:/tools/vcpkg/packages/eigen3_x64-windows/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#target_include_directories(PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(TrafficAssignmentApproaches
    "src/main.cpp"
    # Add other source files as needed
)


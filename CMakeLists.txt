cmake_minimum_required(VERSION 3.10)

if(CMAKE_HOST_WIN32)
    set(CMAKE_TOOLCHAIN_FILE "C:/tools/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
endif()

project(TrafficAssignmentApproaches)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(TrafficAssignmentApproaches
    "src/main.cpp"
)

target_include_directories(TrafficAssignmentApproaches PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/traffic_assignment
    "C:/tools/vcpkg/packages/eigen3_x64-windows/include/eigen3"
)
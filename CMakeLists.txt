cmake_minimum_required(VERSION 3.10)

project(TrafficAssignmentApproaches LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PkgConfig REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS multiprecision)
find_package(NLopt CONFIG REQUIRED)

add_executable(main src/main.cpp)

target_link_libraries(main PRIVATE
    Eigen3::Eigen
    Boost::multiprecision
    NLopt::nlopt
)
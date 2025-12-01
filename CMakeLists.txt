cmake_minimum_required(VERSION 3.10)

project(TrafficAssignmentApproaches)

set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")

find_package(PkgConfig REQUIRED)
pkg_check_modules(NLOPT REQUIRED nlopt)

add_executable(main src/main.cpp)

# Link NLopt to your target
target_include_directories(main PRIVATE ${NLOPT_INCLUDE_DIRS})
target_link_libraries(main PRIVATE ${NLOPT_LIBRARIES})
target_compile_options(main PRIVATE ${NLOPT_CFLAGS_OTHER})

#target_include_directories(TrafficAssignmentApproaches PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}/include
#    ${CMAKE_CURRENT_SOURCE_DIR}/include/traffic_assignment
#    "C:/tools/vcpkg/packages/eigen3_x64-windows/include/eigen3"
#)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
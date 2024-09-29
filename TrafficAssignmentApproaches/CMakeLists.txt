# CMakeList.txt : CMake project for TrafficAssignmentApproaches, include source and define
# project specific logic here.
#

# Убедитесь, что вы подключили vcpkg через CMake
#if(CMAKE_HOST_WIN32)
#    set(CMAKE_TOOLCHAIN_FILE "C:\tools\vcpkg\scripts\buildsystems\vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
#endif()

set(EIGEN3_INCLUDE_DIR "C:/tools/vcpkg/packages/eigen3_x64-windows/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})

#find_package(Eigen3 REQUIRED)
#target_link_libraries(TrafficAssignmentApproaches Eigen3::Eigen)

# Add source to this project's executable.
add_executable (TrafficAssignmentApproaches
    "src/main.cpp"
    "src/DataProcessor.cpp")

target_include_directories(TrafficAssignmentApproaches 
    PRIVATE 
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    )


if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET TrafficAssignmentApproaches PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.

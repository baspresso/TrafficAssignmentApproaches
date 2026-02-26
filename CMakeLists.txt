cmake_minimum_required(VERSION 3.10)

project(TrafficAssignmentApproaches LANGUAGES CXX)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS multiprecision)
find_package(NLopt CONFIG REQUIRED)

add_executable(main src/main.cpp)

target_link_libraries(main PRIVATE
    Eigen3::Eigen
    Boost::multiprecision
    NLopt::nlopt
)

if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_filename_component(MINGW_BIN_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
    foreach(MINGW_RUNTIME_DLL libstdc++-6.dll libgcc_s_seh-1.dll libwinpthread-1.dll)
        if(EXISTS "${MINGW_BIN_DIR}/${MINGW_RUNTIME_DLL}")
            add_custom_command(
                TARGET main
                POST_BUILD
                COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                        "${MINGW_BIN_DIR}/${MINGW_RUNTIME_DLL}"
                        "$<TARGET_FILE_DIR:main>"
            )
        endif()
    endforeach()
endif()

cmake_minimum_required(VERSION 3.14)

project(TrafficAssignmentApproaches LANGUAGES CXX)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

find_package(Eigen3 REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS multiprecision)
find_package(NLopt CONFIG REQUIRED)

# --- OptimLib (header-only-like usage via compiled static library) ---
FetchContent_Declare(
    optimlib
    GIT_REPOSITORY https://github.com/kthohr/optim.git
    GIT_TAG        master
    GIT_SUBMODULES include/BaseMatrixOps
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(optimlib)

add_library(optimlib STATIC
    ${optimlib_SOURCE_DIR}/src/unconstrained/de.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/de_prmm.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/pso.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/pso_dv.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/nm.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/gd.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/bfgs.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/lbfgs.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/cg.cpp
    ${optimlib_SOURCE_DIR}/src/unconstrained/newton.cpp
    ${optimlib_SOURCE_DIR}/src/constrained/sumt.cpp
    ${optimlib_SOURCE_DIR}/src/line_search/more_thuente.cpp
    ${optimlib_SOURCE_DIR}/src/zeros/broyden.cpp
    ${optimlib_SOURCE_DIR}/src/zeros/broyden_df.cpp
)
target_include_directories(optimlib PUBLIC ${optimlib_SOURCE_DIR}/include)
target_compile_definitions(optimlib PUBLIC OPTIM_ENABLE_EIGEN_WRAPPERS)
target_link_libraries(optimlib PUBLIC Eigen3::Eigen)

add_executable(cndp_solver src/cndp_solver.cpp)

target_link_libraries(cndp_solver PRIVATE
    Eigen3::Eigen
    Boost::multiprecision
    NLopt::nlopt
    optimlib
)

add_executable(tap_solver src/tap_solver.cpp)

target_link_libraries(tap_solver PRIVATE
    Eigen3::Eigen
    Boost::multiprecision
)

# MinGW debug builds hit the COFF section limit on large TUs; -Wa,-mbig-obj raises it.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(cndp_solver PRIVATE $<$<CONFIG:Debug>:-Wa,-mbig-obj>)
    target_compile_options(tap_solver PRIVATE $<$<CONFIG:Debug>:-Wa,-mbig-obj>)
endif()

if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_filename_component(MINGW_BIN_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
    foreach(MINGW_RUNTIME_DLL libstdc++-6.dll libgcc_s_seh-1.dll libwinpthread-1.dll)
        if(EXISTS "${MINGW_BIN_DIR}/${MINGW_RUNTIME_DLL}")
            foreach(TARGET_EXE cndp_solver tap_solver)
                add_custom_command(
                    TARGET ${TARGET_EXE}
                    POST_BUILD
                    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                            "${MINGW_BIN_DIR}/${MINGW_RUNTIME_DLL}"
                            "$<TARGET_FILE_DIR:${TARGET_EXE}>"
                )
            endforeach()
        endif()
    endforeach()
endif()

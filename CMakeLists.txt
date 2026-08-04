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
# Boost.Multiprecision is header-only; the Debian Boost package ships only the
# master BoostConfig.cmake (no per-component configs), so request headers only.
find_package(Boost CONFIG REQUIRED)
find_package(NLopt CONFIG REQUIRED)
# Debian's libnlopt-cxx-dev exports NLopt::nlopt_cxx; libnlopt-dev exports
# NLopt::nlopt. Alias to a single name so the target_link_libraries calls
# below don't have to care which one find_package picked up.
if(TARGET NLopt::nlopt_cxx AND NOT TARGET NLopt::nlopt)
    add_library(NLopt::nlopt INTERFACE IMPORTED)
    target_link_libraries(NLopt::nlopt INTERFACE NLopt::nlopt_cxx)
endif()

# --- toml++ (header-only TOML parser) ---
FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG        v3.4.0
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(tomlplusplus)

# --- OptimLib (header-only-like usage via compiled static library) ---
# Pinned to a specific commit for reproducible builds (was: master).
# GIT_SHALLOW is intentionally omitted: a shallow fetch of an arbitrary commit
# SHA is rejected by GitHub (only branch tips are advertised), so a full clone
# is required to check out a pinned hash.
FetchContent_Declare(
    optimlib
    GIT_REPOSITORY https://github.com/kthohr/optim.git
    GIT_TAG        5453f48aca695e6fef123677f35cb4d22e356e73
    GIT_SUBMODULES include/BaseMatrixOps
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
# SYSTEM so -Wall/-Wextra on first-party targets don't flag optim's own headers.
target_include_directories(optimlib SYSTEM PUBLIC ${optimlib_SOURCE_DIR}/include)
target_compile_definitions(optimlib PUBLIC OPTIM_ENABLE_EIGEN_WRAPPERS)
target_link_libraries(optimlib PUBLIC Eigen3::Eigen)
# Required to link the static lib into the Python shared module.
set_target_properties(optimlib PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(cndp_solver src/cndp_solver.cpp)

target_link_libraries(cndp_solver PRIVATE
    Eigen3::Eigen
    Boost::headers
    NLopt::nlopt
    optimlib
    tomlplusplus::tomlplusplus
)

add_executable(tap_solver src/tap_solver.cpp)

target_link_libraries(tap_solver PRIVATE
    Eigen3::Eigen
    Boost::headers
    tomlplusplus::tomlplusplus
)

# Warning flags on first-party targets only (not third-party FetchContent libs).
set(TA_WARNING_FLAGS $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra>)
target_compile_options(cndp_solver PRIVATE ${TA_WARNING_FLAGS})
target_compile_options(tap_solver PRIVATE ${TA_WARNING_FLAGS})

# --- Python bindings (traffic_assignment._core) ---
# Built by scikit-build-core via pyproject.toml; OFF for plain CLI builds.
option(BUILD_PYTHON_BINDINGS "Build the traffic_assignment Python extension module" OFF)

if(BUILD_PYTHON_BINDINGS)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_core bindings/py_bindings.cpp)
    target_link_libraries(_core PRIVATE
        Eigen3::Eigen
        Boost::headers
        NLopt::nlopt
        optimlib
    )
    target_compile_options(_core PRIVATE ${TA_WARNING_FLAGS})
    install(TARGETS _core DESTINATION traffic_assignment COMPONENT python)
endif()

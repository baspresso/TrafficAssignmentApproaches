add_library(traffic_assignment_core STATIC
  src/network.cpp
  src/assignment.cpp
  src/design.cpp
)
add_library(traffic_assignment::core ALIAS traffic_assignment_core)
target_compile_features(traffic_assignment_core PUBLIC cxx_std_20)
set_target_properties(traffic_assignment_core PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_EXTENSIONS OFF
)
target_include_directories(traffic_assignment_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(traffic_assignment_core PRIVATE
  Eigen3::Eigen Boost::headers NLopt::nlopt traffic_assignment_optimlib Threads::Threads
)
set(TA_WARNING_FLAGS $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra>)
target_compile_options(traffic_assignment_core PRIVATE ${TA_WARNING_FLAGS})

if(BUILD_CPP_APPS)
  foreach(app IN ITEMS tap_solver cndp_solver)
    add_executable(${app} apps/${app}.cpp)
    target_link_libraries(${app} PRIVATE traffic_assignment::core tomlplusplus::tomlplusplus)
    target_compile_options(${app} PRIVATE ${TA_WARNING_FLAGS})
    set_target_properties(${app} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
  endforeach()
endif()

if(BUILD_PYTHON_BINDINGS)
  set(PYBIND11_FINDPYTHON ON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core
    bindings/module.cpp bindings/network.cpp bindings/options.cpp
    bindings/results.cpp bindings/solvers.cpp
  )
  target_link_libraries(_core PRIVATE traffic_assignment::core)
  target_compile_options(_core PRIVATE ${TA_WARNING_FLAGS})
  install(TARGETS _core LIBRARY DESTINATION traffic_assignment COMPONENT python
                       RUNTIME DESTINATION traffic_assignment COMPONENT python)
endif()

if(BUILD_NATIVE_TESTS)
  add_executable(test_native_api ../tests/cpp/test_native_api.cpp)
  target_link_libraries(test_native_api PRIVATE traffic_assignment::core)
  add_test(NAME native_api COMMAND test_native_api)

  add_executable(test_route_statistics ../tests/cpp/test_route_statistics.cpp)
  target_include_directories(test_route_statistics PRIVATE src)
  target_link_libraries(test_route_statistics PRIVATE traffic_assignment::core Eigen3::Eigen)
  add_test(NAME route_statistics COMMAND test_route_statistics)
endif()

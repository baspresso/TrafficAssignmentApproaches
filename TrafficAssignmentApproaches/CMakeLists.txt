# CMakeList.txt : CMake project for TrafficAssignmentApproaches, include source and define
# project specific logic here.
#


# Add source to this project's executable.
add_executable (TrafficAssignmentApproaches
    "src/main.cpp"
    )

target_include_directories(TrafficAssignmentApproaches 
    PRIVATE 
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    )
#find_package(Boost REQUIRED)
#target_link_libraries(TrafficAssignmentApproaches Boost::Boost)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET TrafficAssignmentApproaches PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.

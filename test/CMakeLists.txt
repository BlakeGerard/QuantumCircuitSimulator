add_executable(test test.cpp)
target_link_libraries(test qcs_core)
target_link_libraries(test Eigen3::Eigen)
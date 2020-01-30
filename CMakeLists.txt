cmake_minimum_required (VERSION 3.10)
project (qcs)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${qcs_SOURCE_DIR}/include)
include_directories(lib/eigen-3.3.7)

add_subdirectory(src)
add_subdirectory(test)

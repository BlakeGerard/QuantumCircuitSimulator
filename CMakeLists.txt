# Configuration
cmake_minimum_required (VERSION 3.10)
set(CMAKE_CXX_FLAGS "-Wall -g -O3 -fopenmp")
set(CMAKE_CXX_STANDARD 11)
include(GNUInstallDirs)

project(qcs VERSION 1.0)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Catch2 REQUIRED)

add_library(qcs STATIC  
            src/qubit.cpp
            src/gate_matrices.cpp
            src/q_circuit_gates.cpp
            src/q_circuit_ops.cpp)

add_executable(alltests
               src/qubit.cpp 
               src/gate_matrices.cpp
               src/q_circuit_gates.cpp
               src/q_circuit_ops.cpp
               test/test.cpp)

target_include_directories(qcs PRIVATE include)
target_include_directories(alltests PRIVATE include)

target_link_libraries(qcs Eigen3::Eigen)
target_link_libraries(qcs Catch2::Catch2)
target_link_libraries(alltests Eigen3::Eigen)
target_link_libraries(alltests Catch2::Catch2)

set_target_properties(qcs PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(qcs PROPERTIES PUBLIC_HEADER include/qcs.h)

install(TARGETS qcs
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
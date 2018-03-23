cmake_minimum_required(VERSION 3.9)
project(Case1)

set(CMAKE_CXX_STANDARD 11)

add_executable(Case1 main.cpp vtkwrite.cpp callfunctions.h zeros_int.cpp zeros_double.cpp)

find_package(BLAS)
target_link_libraries(Case1 ${BLAS_LIBRARIES})


cmake_minimum_required(VERSION 3.18)
project(blis_eigen_issue C)

find_package(LAPACK REQUIRED)

add_executable(eigen_example main.c)
target_link_libraries(eigen_example PRIVATE LAPACK::LAPACK)

enable_testing()
add_test(NAME eigen_example COMMAND eigen_example)

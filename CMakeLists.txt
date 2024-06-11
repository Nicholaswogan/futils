cmake_minimum_required(VERSION 3.14)

project(FUTILS LANGUAGES Fortran VERSION "0.1.10")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

add_subdirectory(src)
add_subdirectory(tests)



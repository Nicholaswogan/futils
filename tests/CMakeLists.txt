add_executable(test_futils test_futils.f90)
target_link_libraries(test_futils futils)
target_include_directories(test_futils PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

# Special compiler options
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_futils PRIVATE -Wunused -fcheck=all)
  endif()
endif()
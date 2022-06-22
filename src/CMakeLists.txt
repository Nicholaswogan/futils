add_library(futils_fastmath futils_fastmath.f90)

add_library(futils futils.f90 futils_brent.f90)
target_link_libraries(futils futils_fastmath)

# Special compiler options
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(futils_fastmath PRIVATE -Wunused -fcheck=all)
    target_compile_options(futils PRIVATE -Wunused -fcheck=all)
  endif()
  
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(futils_fastmath PRIVATE -ffast-math -finline-functions -ftree-vectorize -funroll-loops -march=native)
  endif()
endif()


add_library(futils
  futils_brent.f90 
  futils_rebin.f90
  futils_special.f90
  futils.f90 
)

# Special compiler options
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(futils PRIVATE -Wall -fcheck=all)
  endif()
endif()

add_library(specfun specfun.f)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(specfun PRIVATE -std=legacy)
endif()

add_library(futils
  futils_brent.f90 
  futils_rebin.f90
  futils_special.f90
  futils_misc.f90
  futils.f90 
)
target_link_libraries(futils specfun)

# Special compiler options
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(futils PRIVATE -Wunused -Wimplicit-interface -fimplicit-none -fcheck=all,no-array-temps)
  endif()
endif()

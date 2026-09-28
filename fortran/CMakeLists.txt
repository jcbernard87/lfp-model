# fortran/lfp.f90 is added in T5.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/lfp.f90)
  add_executable(lfp_f lfp.f90)
  target_link_libraries(lfp_f PRIVATE bandsolver_f)
  target_include_directories(lfp_f PRIVATE ${CMAKE_BINARY_DIR}/modules)
  # No fused multiply-add contraction: faithful mode reproduces the original program exactly
  target_compile_options(lfp_f PRIVATE $<$<Fortran_COMPILER_ID:GNU>:-ffp-contract=off -std=f2008 -Wall>)
endif()

# cpp/lfp.cpp is a single self-contained file with its own embedded BAND solver.
add_executable(lfp_cpp lfp.cpp)
# No fused multiply-add contraction: faithful mode reproduces the original program exactly
target_compile_options(lfp_cpp PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-ffp-contract=off -Wall -Wextra>)

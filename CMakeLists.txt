cmake_minimum_required(VERSION 3.20)
project(lfp_model VERSION 0.3.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The Fortran program uses bandsolver. By default it is fetched from GitHub;
# use a local checkout with -DFETCHCONTENT_SOURCE_DIR_BANDSOLVER=/path/to/bandsolver
include(FetchContent)
FetchContent_Declare(bandsolver
  GIT_REPOSITORY https://github.com/jcbernard87/bandsolver.git
  GIT_TAG v0.1.1)
set(BANDSOLVER_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BANDSOLVER_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(bandsolver)

enable_testing()
add_subdirectory(fortran)
add_subdirectory(cpp)

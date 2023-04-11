# Create tests for header-only files in namespace hd

##################################################################
# usage: execute "cmake .. && make" from build directory
#    or: execute "cmake .. -G"Ninja" && ninja from build directory
#
# alternative to build:  "cmake --build ." from build directory
##################################################################
#]]
#[[ hint: backet comment will not be formated by cmake-format]]

cmake_minimum_required(VERSION 3.25)

# to use gcc-12 and g++-12 instead of default c / c++ (must come before project
# line [[
# set(CMAKE_C_COMPILER "gcc-12")
# set(CMAKE_CXX_COMPILER "g++-12")
# ]]

# CMAKE_TOOLCHAIN_FILE must be defined BEFORE project(...)!
if(WIN32)
  set(CMAKE_TOOLCHAIN_FILE "C:/Users/danie/prg/cpp/vcpkg/scripts/buildsystems/vcpkg.cmake")
else()
  set(CMAKE_TOOLCHAIN_FILE "/Users/daniel/prg/cpp/vcpkg/scripts/buildsystems/vcpkg.cmake")
endif()

#
# for trial of new features, libs or system values"
#
project(
  experiments
  VERSION 0.1
  LANGUAGES CXX)

# set(CMAKE_VERBOSE_MAKEFILE ON)  // alternative: make VERBOSE=1
# add_definitions(-Dxxx)

# make sure cmake is called with indication of CMAKE_TOOLCHAIN_FILE (vcpkg)
# e.g. via .zshrc export CMAKE_TOOLCHAIN_FILE="/Users/daniel/prg/cpp/vcpkg/scripts/buildsystems/vcpkg.cmake"
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  message(FATAL_ERROR
          "CMAKE_TOOLCHAIN_FILE must be defined to use vcpkg. Configuration aborted.")
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# set a default build type: Debug | RelWithDebInfo | Release | MinSizeRel
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
  message(STATUS "Build type not specified: use Debug mode by default.")
endif()

if(CMAKE_BUILD_TYPE EQUAL "Debug")
  if(MSVC)
    # warning level 4 and all warnings as errors
    add_compile_options(/W4 /WX)
  else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic -Werror)
  endif()
endif()

find_package(doctest REQUIRED)

# incrementally add test files needed
add_executable(hd_functions_test hd_functions_test.cpp)     #dep: ...

target_link_libraries(hd_functions_test PRIVATE doctest::doctest)
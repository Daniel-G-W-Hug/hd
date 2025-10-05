# Create tests for header-only files in namespace hd

##################################################################
# usage: execute "cmake .. && make" from build directory
#    or: execute "cmake .. -G"Ninja" && ninja from build directory
#
# alternative to build:  "cmake --build ." from build directory
##################################################################
#]]
#[[ hint: backet comment will not be formated by cmake-format]]

cmake_minimum_required(VERSION 3.29)

#
# for trial of new features, libs or system values"
#
project(
  hd_lib
  VERSION 0.1
  LANGUAGES CXX)

# set(CMAKE_VERBOSE_MAKEFILE ON)  // alternative: make VERBOSE=1
# add_definitions(-Dxxx)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# set a default build type: Debug | RelWithDebInfo | Release | MinSizeRel
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
  message(STATUS "Build type not specified: use Debug mode by default.")
endif()

if(MSVC)
  # warning level 4 and all warnings as errors
  add_compile_options(/W4 /WX)
else()
  # lots of warnings and all warnings as errors
  #add_compile_options(-Wall -Wextra -pedantic -Werror)
  add_compile_options(-Wall -Wextra -pedantic)
endif()

find_package(doctest REQUIRED)
#find_package(date REQUIRED)
find_package(fmt CONFIG REQUIRED)

add_executable(hd_functions_test hd_functions_test.cpp)   #dep: doctest
target_link_libraries(hd_functions_test PRIVATE doctest::doctest)

add_executable(hd_solver_test hd_solver_test.cpp)   #dep: doctest
target_link_libraries(hd_solver_test PRIVATE doctest::doctest)
target_include_directories(hd_solver_test PRIVATE ${CMAKE_SOURCE_DIR}/../../include/mdspan/include)
target_compile_definitions(hd_solver_test PUBLIC MDSPAN_USE_BRACKET_OPERATOR=1 MDSPAN_CXX_STANDARD=23)
target_link_libraries(hd_solver_test PRIVATE fmt::fmt-header-only)

add_executable(hd_determinant_test hd_determinant_test.cpp)   #dep: doctest
target_link_libraries(hd_determinant_test PRIVATE doctest::doctest)
target_include_directories(hd_determinant_test PRIVATE ${CMAKE_SOURCE_DIR}/../../include/mdspan/include)
target_compile_definitions(hd_determinant_test PUBLIC MDSPAN_USE_BRACKET_OPERATOR=1 MDSPAN_CXX_STANDARD=23)
target_link_libraries(hd_determinant_test PRIVATE fmt::fmt-header-only)

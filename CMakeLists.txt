# Create tests for header-only files in namespace hd

##################################################################
# usage: execute "cmake .. && make" from build directory
#    or: execute "cmake .. -G"Ninja" && ninja from build directory
#
# alternative to build:  "cmake --build ." from build directory
#
# vcpkg usage on Windows:
#    cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
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

# Cross-platform dependency management with vcpkg auto-detection

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

# Dependency management - detect vcpkg and use appropriate method
if(WIN32 AND DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg")
  message(STATUS "Detected vcpkg toolchain on Windows")

  # Using vcpkg - CONFIG mode find_package calls
  find_package(doctest CONFIG REQUIRED)
  if(doctest_FOUND)
    message(STATUS "✓ Found vcpkg doctest")
  else()
    message(FATAL_ERROR "doctest not found via vcpkg - run: vcpkg install doctest")
  endif()

  find_package(fmt CONFIG REQUIRED)
  if(fmt_FOUND)
    message(STATUS "✓ Found vcpkg fmt: ${fmt_VERSION}")
  else()
    message(FATAL_ERROR "fmt not found via vcpkg - run: vcpkg install fmt")
  endif()

  #find_package(date CONFIG REQUIRED)  # Uncomment if needed

  # For mdspan with vcpkg, try to find it, fallback to local copy
  find_package(mdspan CONFIG QUIET)
  if(NOT mdspan_FOUND)
    message(STATUS "mdspan not found in vcpkg, using local copy at ../../include/mdspan/include")
    set(MDSPAN_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../include/mdspan/include")
  else()
    message(STATUS "✓ Found vcpkg mdspan")
  endif()
else()
  # Traditional approach (macOS with brew, or manual installation)
  message(STATUS "Using traditional dependency management")
  find_package(doctest REQUIRED)
  find_package(fmt CONFIG REQUIRED)
  #find_package(date REQUIRED)  # Uncomment if needed

  # Always use local mdspan copy for non-vcpkg builds
  set(MDSPAN_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../include/mdspan/include")
  message(STATUS "Using local mdspan at ${MDSPAN_INCLUDE_DIR}")
endif()

# Function to configure test executable with proper dependencies
function(configure_test_executable target_name source_file needs_mdspan needs_fmt)
  add_executable(${target_name} ${source_file})

  # Link doctest (same interface for both vcpkg and traditional)
  target_link_libraries(${target_name} PRIVATE doctest::doctest)

  # Configure mdspan if needed
  if(needs_mdspan)
    if(WIN32 AND DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg" AND mdspan_FOUND)
      target_link_libraries(${target_name} PRIVATE std::mdspan)
    else()
      target_include_directories(${target_name} PRIVATE ${MDSPAN_INCLUDE_DIR})
    endif()
    target_compile_definitions(${target_name} PUBLIC MDSPAN_USE_BRACKET_OPERATOR=1 MDSPAN_CXX_STANDARD=23)
  endif()

  # Configure fmt if needed
  if(needs_fmt)
    if(WIN32 AND DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg")
      target_link_libraries(${target_name} PRIVATE fmt::fmt)
    else()
      target_link_libraries(${target_name} PRIVATE fmt::fmt-header-only)
    endif()
  endif()
endfunction()

# Configure test executables
configure_test_executable(hd_functions_test hd_functions_test.cpp FALSE FALSE)
configure_test_executable(hd_solver_test hd_solver_test.cpp TRUE TRUE)
configure_test_executable(hd_determinant_test hd_determinant_test.cpp TRUE TRUE)

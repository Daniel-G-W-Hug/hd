

# incrementally add test files needed

add_executable(hd_ga_test hd_ga_test.cpp)                 #dep: doctest, fmt

find_package(doctest REQUIRED)
find_package(fmt CONFIG REQUIRED)

target_link_libraries(hd_ga_test PRIVATE doctest::doctest fmt::fmt-header-only)


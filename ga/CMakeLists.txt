# incrementally add test files needed

set(EXEC_NAME1 hd_ga_test)
set(EXEC_NAME2 hd_ga_lua)

add_executable(${EXEC_NAME1} hd_ga_test.cpp)                 #dep: doctest, fmt
add_executable(${EXEC_NAME2} hd_ga_lua.cpp)                  #dep: doctest, fmt, lua, sol2

target_include_directories(${EXEC_NAME2} PRIVATE ${CMAKE_SOURCE_DIR}/../../include/sol2/include)

find_package(doctest REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(Lua 5.4 REQUIRED)

target_link_libraries(${EXEC_NAME1} PRIVATE doctest::doctest fmt::fmt-header-only)

target_include_directories(${EXEC_NAME2} PRIVATE ${LUA_INCLUDE_DIR})
target_link_libraries(${EXEC_NAME2} PRIVATE doctest::doctest fmt::fmt-header-only ${LUA_LIBRARIES})

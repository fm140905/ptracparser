cmake_minimum_required(VERSION 3.16.3)
project(PTRACParser CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

enable_testing()
add_subdirectory(googletest)
add_subdirectory(src)
add_subdirectory(tests)

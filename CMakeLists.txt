cmake_minimum_required(VERSION 3.16.3)
project(PTRACParser CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

option(ENABLE_UNIT_TESTS "Enable unit tests" OFF)
message(STATUS "Enable testing: ${ENABLE_UNIT_TESTS}")

# Import FetchContent module:
include(FetchContent)
if (ENABLE_UNIT_TESTS)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.11.0
    )
    
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    add_subdirectory(tests)
endif()

FetchContent_Declare(
    highfive
    URL
    https://github.com/BlueBrain/HighFive/archive/refs/tags/v2.7.1.zip
)
FetchContent_MakeAvailable(highfive)

add_subdirectory(src)

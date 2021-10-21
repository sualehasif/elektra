cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

project(
  elektra
  VERSION 0.1.0
  DESCRIPTION "A collection of parallel dynamic graph algorithms"
  LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set up the CMAKE Generator
set(CMAKE_GENERATOR "Ninja")
set(CMAKE_MAKE_PROGRAM "Ninja")

include(CheckCXXCompilerFlag)
include(GNUInstallDirs)

# Setting up the C++ CMAKE_CXX_FLAGS initial
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")

# Setting up some unused checks
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} \
    -Wno-unused-parameter -Wno-unused-variable \
    -Wno-unused-function -Wno-unused-private-field \
    -Wno-unused-local-typedef -Wno-unused-const-variable \
    -Wno-unused-value \
    ")

# # Compile with no parallelism if debugging or asked to be linear
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARLAY_SEQUENTIAL")
endif()

message(STATUS "--------------- General configuration -------------")
message(STATUS "CMake Generator: ${CMAKE_GENERATOR}")
message(STATUS "Compiler: \
    ${CMAKE_CXX_COMPILER_ID} \
    ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Build type: \
    ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS:
    ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG:
    ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE:
    ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO:
    ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS
    ${CMAKE_CXX_LINKER_FLAGS}")
message(STATUS "CMAKE_INSTALL_PREFIX: \
    ${CMAKE_INSTALL_PREFIX}")

include(CTest)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/elektra)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/parlaylib/include/)
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googlemock/include)
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

link_directories(${Boost_LIBRARY_DIRS})

enable_testing()

add_executable(elektra elektra.cpp)

target_link_libraries(elektra ${Boost_LIBRARIES})

# Testing
add_executable(elektra_test test/elektra_test.cpp)
target_link_libraries(elektra_test gmock_main gtest_main)
include(GoogleTest)
gtest_discover_tests(elektra_test)

# Set up the CMAKE Generator
set(CMAKE_GENERATOR "Ninja")
set(CMAKE_MAKE_PROGRAM "Ninja")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

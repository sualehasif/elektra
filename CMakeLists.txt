cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

project(elektra VERSION 0.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(CTest)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/elektra)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/parlaylib/include/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/include)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/googletest/)

find_package(Boost REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

link_directories(${Boost_LIBRARY_DIRS})

message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")

enable_testing()

add_executable(elektra elektra.cpp)

TARGET_LINK_LIBRARIES(elektra ${Boost_LIBRARIES})

# Testing
add_executable(elektra_test test/elektra_test.cpp)
target_link_libraries(
  elektra_test
  gtest_main
)
include(GoogleTest)
gtest_discover_tests(elektra_test)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

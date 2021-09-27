cmake_minimum_required(VERSION 3.0.0)
project(elektra VERSION 0.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(CTest)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/elektra)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/parlaylib/include/)

find_package(Boost REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")

enable_testing()

add_executable(elektra elektra.cpp)

TARGET_LINK_LIBRARIES(elektra ${Boost_LIBRARIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

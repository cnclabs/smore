# file: CMakeLists.txt
# date: 2022-10-07


cmake_minimum_required(VERSION 3.14)
set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

project(smore LANGUAGES C CXX)


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  error("Argument '-fopenmp' only supported by GCC on mac!")
else()
  set(CMAKE_CXX_FLAGS "-pthread -fopenmp -std=c++14 -funroll-loops -O3 -march=native")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")


include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/models)

file(GLOB CPP_SOURCES
  "${PROJECT_SOURCE_DIR}/src/*.cpp" 
  "${PROJECT_SOURCE_DIR}/src/model/*.cpp"
)

add_library(smore STATIC ${CPP_SOURCES})


add_subdirectory(${PROJECT_SOURCE_DIR}/cli)

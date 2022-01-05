cmake_minimum_required(VERSION 3.10)
project(
  baobzi
  LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

set(BAOBZI_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/extern/msgpack-c/include
  ${EIGEN3_INCLUDE_DIR}
  )


file(GLOB LIB_SOURCES "src/*.cpp")
add_library(baobzi SHARED ${LIB_SOURCES})
target_include_directories(baobzi PUBLIC ${BAOBZI_INCLUDES})
target_link_libraries(baobzi PUBLIC)

set(EXAMPLE_SOURCE_CPP "examples/baobzi_timing.cpp")
add_executable(baobzi_timing_cpp ${EXAMPLE_SOURCE_CPP})
target_include_directories(baobzi_timing_cpp PUBLIC ${BAOBZI_INCLUDES})
target_link_libraries(baobzi_timing_cpp PUBLIC OpenMP::OpenMP_CXX)

set(EXAMPLE_SOURCE_C "examples/baobzi_timing.c")
add_executable(baobzi_timing_c ${EXAMPLE_SOURCE_C})
target_include_directories(baobzi_timing_c PUBLIC ${BAOBZI_INCLUDES})
target_link_libraries(baobzi_timing_c PUBLIC baobzi OpenMP::OpenMP_CXX m)

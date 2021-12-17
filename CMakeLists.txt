cmake_minimum_required(VERSION 3.10)
project(
  baobzi
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

set(BAOBZI_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/extern/msgpack-c/include
  ${EIGEN3_INCLUDE_DIR}
  )

file(GLOB EXAMPLE_SOURCES "examples/baobzi_timing.cpp")

add_executable(baobzi_timing ${EXAMPLE_SOURCES})
target_include_directories(baobzi_timing PUBLIC ${BAOBZI_INCLUDES})
target_link_libraries(baobzi_timing PUBLIC OpenMP::OpenMP_CXX)

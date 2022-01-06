cmake_minimum_required(VERSION 3.10)
project(
  baobzi
  LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

set(BAOBZI_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/extern/msgpack-c/include
  ${PROJECT_SOURCE_DIR}/extern/eigen
  )

option(BUILD_STATIC "Build the static library" OFF)
option(BUILD_SHARED "Build the shared library" ON)
option(BUILD_EXAMPLES "Build C/C++ examples" ON)

file(GLOB LIB_SOURCES_GENERIC "src/*0.cpp")
add_library(baobzi_generic OBJECT ${LIB_SOURCES_GENERIC})
target_include_directories(baobzi_generic PRIVATE ${BAOBZI_INCLUDES})
target_compile_options(baobzi_generic PRIVATE)
set_property(TARGET baobzi_generic PROPERTY POSITION_INDEPENDENT_CODE ON)

file(GLOB LIB_SOURCES_AVX "src/*1.cpp")
add_library(baobzi_avx OBJECT ${LIB_SOURCES_AVX})
target_include_directories(baobzi_avx PRIVATE ${BAOBZI_INCLUDES})
target_compile_options(baobzi_avx PRIVATE -mavx)
set_property(TARGET baobzi_avx PROPERTY POSITION_INDEPENDENT_CODE ON)

file(GLOB LIB_SOURCES_AVX2 "src/*2.cpp")
add_library(baobzi_avx2 OBJECT ${LIB_SOURCES_AVX2})
target_include_directories(baobzi_avx2 PRIVATE ${BAOBZI_INCLUDES})
target_compile_options(baobzi_avx2 PRIVATE -mavx2 -mfma)
set_property(TARGET baobzi_avx2 PROPERTY POSITION_INDEPENDENT_CODE ON)

file(GLOB LIB_SOURCES_AVX512 "src/*3.cpp")
add_library(baobzi_avx512 OBJECT ${LIB_SOURCES_AVX512})
target_include_directories(baobzi_avx512 PRIVATE ${BAOBZI_INCLUDES})
target_compile_options(baobzi_avx512 PRIVATE -mavx512f -mfma)
set_property(TARGET baobzi_avx512 PROPERTY POSITION_INDEPENDENT_CODE ON)

if (BUILD_SHARED)
  add_library(baobzi SHARED src/baobzi.cpp)
  target_include_directories(baobzi PUBLIC ${BAOBZI_INCLUDES})
  target_link_libraries(baobzi PUBLIC baobzi_generic baobzi_avx baobzi_avx2 baobzi_avx512)
endif()

if (BUILD_STATIC)
  add_library(baobzi_static STATIC)
  target_link_libraries(baobzi_static PUBLIC baobzi_generic baobzi_avx baobzi_avx2 baobzi_avx512)
  set_target_properties(baobzi_static PROPERTIES OUTPUT_NAME baobzi)
endif()

if (${BUILD_SHARED} AND ${BUILD_EXAMPLES})
  set(EXAMPLE_SOURCE_CPP "examples/baobzi_timing.cpp")
  add_executable(baobzi_timing_cpp ${EXAMPLE_SOURCE_CPP})
  target_include_directories(baobzi_timing_cpp PUBLIC ${BAOBZI_INCLUDES})
  target_link_libraries(baobzi_timing_cpp PUBLIC OpenMP::OpenMP_CXX)
  target_compile_options(baobzi_timing_cpp PRIVATE -march=native)

  set(EXAMPLE_SOURCE_C "examples/baobzi_timing.c")
  add_executable(baobzi_timing_c ${EXAMPLE_SOURCE_C})
  target_include_directories(baobzi_timing_c PUBLIC ${BAOBZI_INCLUDES})
  target_link_libraries(baobzi_timing_c PUBLIC baobzi OpenMP::OpenMP_CXX m)
endif()

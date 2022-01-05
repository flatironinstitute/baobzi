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


file(GLOB LIB_SOURCES_GENERIC "src/*0.cpp")
add_library(baobzi_avx OBJECT ${LIB_SOURCES_GENERIC})
target_include_directories(baobzi_avx PRIVATE ${BAOBZI_INCLUDES})
target_compile_options(baobzi_avx PRIVATE -mavx)
set_property(TARGET baobzi_avx PROPERTY POSITION_INDEPENDENT_CODE ON)

file(GLOB LIB_SOURCES_AVX2 "src/*1.cpp")
add_library(baobzi_avx2 OBJECT ${LIB_SOURCES_AVX2})
target_include_directories(baobzi_avx2 PRIVATE ${BAOBZI_INCLUDES})
target_compile_options(baobzi_avx2 PRIVATE -mavx2 -mfma)
set_property(TARGET baobzi_avx2 PROPERTY POSITION_INDEPENDENT_CODE ON)

file(GLOB LIB_SOURCES_AVX512 "src/*2.cpp")
add_library(baobzi_avx512 OBJECT ${LIB_SOURCES_AVX512})
target_include_directories(baobzi_avx512 PRIVATE ${BAOBZI_INCLUDES})
target_compile_options(baobzi_avx512 PRIVATE -mavx512f -mfma)
set_property(TARGET baobzi_avx512 PROPERTY POSITION_INDEPENDENT_CODE ON)


add_library(baobzi SHARED src/baobzi.cpp)
target_include_directories(baobzi PUBLIC ${BAOBZI_INCLUDES})
target_link_libraries(baobzi PUBLIC baobzi_avx baobzi_avx2 baobzi_avx512)

# add_library(baobzi_static STATIC ${LIB_SOURCES})
# target_include_directories(baobzi_static PUBLIC ${BAOBZI_INCLUDES})
# target_link_libraries(baobzi_static PUBLIC)
# set_target_properties(baobzi_static PROPERTIES OUTPUT_NAME baobzi)


set(EXAMPLE_SOURCE_CPP "examples/baobzi_timing.cpp")
add_executable(baobzi_timing_cpp ${EXAMPLE_SOURCE_CPP})
target_include_directories(baobzi_timing_cpp PUBLIC ${BAOBZI_INCLUDES})
target_link_libraries(baobzi_timing_cpp PUBLIC OpenMP::OpenMP_CXX)

set(EXAMPLE_SOURCE_C "examples/baobzi_timing.c")
add_executable(baobzi_timing_c ${EXAMPLE_SOURCE_C})
target_include_directories(baobzi_timing_c PUBLIC ${BAOBZI_INCLUDES})
target_link_libraries(baobzi_timing_c PUBLIC baobzi OpenMP::OpenMP_CXX m)

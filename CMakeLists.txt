cmake_minimum_required(VERSION 3.5)
project(
  baobzi
  LANGUAGES C CXX)
enable_language(Fortran OPTIONAL)

set(CMAKE_CXX_STANDARD 17)

set (default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message (STATUS "No build type specified. Setting build type to Release.")
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Valid options: Debug, RelWithDebInfo, Release" FORCE)
endif()

set(BAOBZI_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/extern/msgpack-c/include
  ${PROJECT_SOURCE_DIR}/extern/eigen
  ${PROJECT_SOURCE_DIR}/extern/catch2/src
  )

option(BAOBZI_BUILD_STATIC "Build the static library" OFF)
option(BAOBZI_BUILD_SHARED "Build the shared library" ON)
option(BAOBZI_BUILD_EXAMPLES "Build C/C++ examples" ON)
option(BAOBZI_BUILD_TESTS "Build tests" ON)
option(BAOBZI_BUILD_MATLAB "Build MATLAB bindings" OFF)
option(BAOBZI_CPU_DISPATCH "Use CPU dispatch for generic binary" ON)


if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "arm64")
  set(BAOBZI_CPU_DISPATCH OFF)
endif()

if (BAOBZI_CPU_DISPATCH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBAOBZI_CPU_DISPATCH")
endif()

if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "arm64" AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
  set(NATIVE_ARCH_FLAG -mcpu=apple-m1)
else()
  set(NATIVE_ARCH_FLAG -march=native)
endif()


file(GLOB LIB_SOURCES_GENERIC "src/*0.cpp")
add_library(baobzi_generic OBJECT ${LIB_SOURCES_GENERIC})
target_include_directories(baobzi_generic PUBLIC ${BAOBZI_INCLUDES})
target_compile_options(baobzi_generic PRIVATE)
if (NOT BAOBZI_CPU_DISPATCH)
  target_compile_options(baobzi_generic PRIVATE ${NATIVE_ARCH_FLAG})
endif()

set_property(TARGET baobzi_generic PROPERTY POSITION_INDEPENDENT_CODE ON)

if (BAOBZI_CPU_DISPATCH)
  file(GLOB LIB_SOURCES_AVX "src/*1.cpp")
  add_library(baobzi_avx OBJECT ${LIB_SOURCES_AVX})
  target_include_directories(baobzi_avx PUBLIC ${BAOBZI_INCLUDES})
  target_compile_options(baobzi_avx PRIVATE -mavx)
  set_property(TARGET baobzi_avx PROPERTY POSITION_INDEPENDENT_CODE ON)

  file(GLOB LIB_SOURCES_AVX2 "src/*2.cpp")
  add_library(baobzi_avx2 OBJECT ${LIB_SOURCES_AVX2})
  target_include_directories(baobzi_avx2 PUBLIC ${BAOBZI_INCLUDES})
  target_compile_options(baobzi_avx2 PRIVATE -mavx2 -mfma)
  set_property(TARGET baobzi_avx2 PROPERTY POSITION_INDEPENDENT_CODE ON)

  file(GLOB LIB_SOURCES_AVX512 "src/*3.cpp")
  add_library(baobzi_avx512 OBJECT ${LIB_SOURCES_AVX512} )
  target_include_directories(baobzi_avx512 PUBLIC ${BAOBZI_INCLUDES})
  target_compile_options(baobzi_avx512 PRIVATE -mavx512f -mfma)
  set_property(TARGET baobzi_avx512 PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

set(BAOBZI_BINDING_SOURCES src/baobzi.cpp)

if (DEFINED CMAKE_Fortran_COMPILER)
  list(APPEND BAOBZI_BINDING_SOURCES src/baobzi.f90)
endif()

if (BAOBZI_BUILD_SHARED)
  add_library(baobzi SHARED ${BAOBZI_BINDING_SOURCES})
  target_include_directories(baobzi PUBLIC ${BAOBZI_INCLUDES})
  if (BAOBZI_CPU_DISPATCH)
    target_link_libraries(baobzi PRIVATE baobzi_generic baobzi_avx baobzi_avx2 baobzi_avx512)
  else()
    target_link_libraries(baobzi PRIVATE baobzi_generic)
  endif()

  list(APPEND INSTALL_TARGETS baobzi)
endif()

if (BAOBZI_BUILD_STATIC)
  add_library(baobzi_static STATIC ${BAOBZI_BINDING_SOURCES})
  if (BAOBZI_CPU_DISPATCH)
    target_link_libraries(baobzi_static PRIVATE baobzi_generic baobzi_avx baobzi_avx2 baobzi_avx512)
  else()
    target_link_libraries(baobzi_static PRIVATE baobzi_generic)
  endif()

  set_target_properties(baobzi_static PROPERTIES OUTPUT_NAME baobzi)
  list(APPEND INSTALL_TARGETS baobzi_static)
endif()

if (${BAOBZI_BUILD_SHARED} AND ${BAOBZI_BUILD_EXAMPLES})
  set(EXAMPLE_SOURCE_CPP "examples/c++/baobzi_timing.cpp")
  add_executable(baobzi_timing_cpp ${EXAMPLE_SOURCE_CPP})
  target_include_directories(baobzi_timing_cpp PUBLIC ${BAOBZI_INCLUDES})
  target_link_libraries(baobzi_timing_cpp PUBLIC)
  target_compile_options(baobzi_timing_cpp PUBLIC ${NATIVE_ARCH_FLAG})

  set(EXAMPLE_SOURCE_C "examples/C/baobzi_timing.c")
  add_executable(baobzi_timing_c ${EXAMPLE_SOURCE_C})
  target_include_directories(baobzi_timing_c PUBLIC ${BAOBZI_INCLUDES})
  target_link_libraries(baobzi_timing_c PUBLIC baobzi ${NATIVE_ARCH_FLAG} m)
endif()

if (BAOBZI_BUILD_TESTS)
  add_subdirectory(extern/catch2)

  add_executable(test_c tests/test_c.cpp)
  target_link_libraries(test_c PUBLIC Catch2::Catch2WithMain baobzi)

  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/extern/catch2/extras)
  include(CTest)
  include(Catch)
  catch_discover_tests(test_c)
endif()

if (BAOBZI_BUILD_MATLAB)
  find_package(Matlab REQUIRED)
  matlab_add_mex(NAME baobzi_mex SRC src/matlab/baobzi_mex.cpp LINK_TO baobzi)
  target_include_directories(baobzi_mex PUBLIC ${BAOBZI_INCLUDES})
  install(FILES ${PROJECT_SOURCE_DIR}/src/matlab/baobzi.m
    DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/matlab)
  install(TARGETS baobzi_mex DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/matlab)
endif()

include(GNUInstallDirs)
install(TARGETS ${INSTALL_TARGETS})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${PROJECT_SOURCE_DIR}/src/julia/baobzi.jl
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/julia)
install(FILES ${PROJECT_SOURCE_DIR}/src/python/baobzi.py
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/python)
install(FILES ${PROJECT_SOURCE_DIR}/src/baobzi.f90
  ${CMAKE_CURRENT_BINARY_DIR}/baobzi.mod
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/fortran)

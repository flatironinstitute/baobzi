cmake_minimum_required(VERSION 3.15)
project(
  baobzi
  LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

set (default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message (STATUS "No build type specified. Setting build type to Release.")
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Valid options: Debug, RelWithDebInfo, Release" FORCE)
endif()

set(BAOBZI_INCLUDES
  ${PROJECT_SOURCE_DIR}/include
)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)
if (ipo_supported)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(WARNING "IPO is not supported: ${error}")
endif()

include(FetchContent)
FetchContent_Declare(
    polyfit
    GIT_REPOSITORY https://github.com/blackwer/polyfit.git
    GIT_TAG        ND-coeff-access
  )
set(Build_Tests OFF CACHE BOOL "" FORCE) # gtest... tests
set(MONOFIT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(MONOFIT_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(polyfit)

option(BAOBZI_BUILD_STATIC "Build the static library" OFF)
option(BAOBZI_BUILD_SHARED "Build the shared library" ON)
option(BAOBZI_BUILD_EXAMPLES "Build C/C++ examples" ON)
option(BAOBZI_BUILD_TESTS "Build tests" ON)
option(BAOBZI_BUILD_MATLAB "Build MATLAB bindings" OFF)
option(BAOBZI_BUILD_FORTRAN "Build Fortran bindings" OFF)
option(BAOBZI_CPU_DISPATCH "Use CPU dispatch for generic binary" OFF)

# Needed for ND polyfit performance
add_compile_options(-funroll-loops -funsafe-math-optimizations)
add_link_options(-funroll-loops -funsafe-math-optimizations)

if (BAOBZI_CPU_DISPATCH)
  add_compile_options(-DBAOBZI_CPU_DISPATCH)

  function (add_dispatch_library isa_name)
    set(name baobzi_${isa_name})
    add_library(${name} SHARED src/baobzi_dispatch.cpp)
    target_include_directories(${name} PUBLIC ${BAOBZI_INCLUDES})
    set_property(TARGET ${name} PROPERTY POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${name} PRIVATE polyfit::polyfit)
    if (NOT MSVC)
      if (${isa_name} STREQUAL "x86_64")
        target_compile_options(${name} PRIVATE -march=x86-64 -DBAOBZI_ISA=0)
      elseif (${isa_name} STREQUAL "x86_64_v2")
        target_compile_options(${name} PRIVATE -march=x86-64-v2 -DBAOBZI_ISA=1)
      elseif(${isa_name} STREQUAL "x86_64_v3")
        target_compile_options(${name} PRIVATE -march=x86-64-v3 -DBAOBZI_ISA=2)
      elseif(${isa_name} STREQUAL "x86_64_v4")
        target_compile_options(${name} PRIVATE -march=x86-64-v4 -DBAOBZI_ISA=3)
      endif()
    else()
      message(FATAL_ERROR "MSVC dispatch not implemented yet")
    endif()
  endfunction()

  add_dispatch_library(x86_64)
  add_dispatch_library(x86_64_v2)
  add_dispatch_library(x86_64_v3)
  add_dispatch_library(x86_64_v4)
endif()

set(BAOBZI_BINDING_SOURCES src/baobzi.cpp)
if (${BAOBZI_BUILD_FORTRAN})
  enable_language(Fortran)
  if (DEFINED CMAKE_Fortran_COMPILER)
    list(APPEND BAOBZI_BINDING_SOURCES src/baobzi.f90)
    install(FILES ${PROJECT_SOURCE_DIR}/src/baobzi.f90
      ${CMAKE_CURRENT_BINARY_DIR}/baobzi.mod
      DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/fortran)
  endif()
endif()

if (${BAOBZI_BUILD_SHARED})
  add_library(baobzi SHARED ${BAOBZI_BINDING_SOURCES})
  target_include_directories(baobzi PUBLIC ${BAOBZI_INCLUDES})
  if (BAOBZI_CPU_DISPATCH)
    target_link_libraries(baobzi PRIVATE polyfit::polyfit dl)
  else()
    target_link_libraries(baobzi PRIVATE polyfit::polyfit)
  endif()

  list(APPEND INSTALL_TARGETS baobzi)
endif()

if (${BAOBZI_BUILD_STATIC})
  add_library(baobzi_static STATIC ${BAOBZI_BINDING_SOURCES})
  if (BAOBZI_CPU_DISPATCH)
    message(FATAL_ERROR "Static library with dispatch not implemented")
  else()
    target_link_libraries(baobzi_static PRIVATE polyfit::polyfit)
  endif()

  set_target_properties(baobzi_static PROPERTIES OUTPUT_NAME baobzi)
  list(APPEND INSTALL_TARGETS baobzi_static)
endif()

if (BAOBZI_BUILD_EXAMPLES)
  add_executable(baobzi_timing_cpp "examples/c++/baobzi_timing.cpp")
  target_include_directories(baobzi_timing_cpp PUBLIC ${BAOBZI_INCLUDES})
  target_link_libraries(baobzi_timing_cpp PUBLIC polyfit::polyfit)

  if (BAOBZI_BUILD_SHARED)
    add_executable(baobzi_timing_c "examples/C/baobzi_timing.c")
    target_include_directories(baobzi_timing_c PUBLIC ${BAOBZI_INCLUDES})
    target_link_libraries(baobzi_timing_c PUBLIC baobzi m)
  endif()
endif()

if (BAOBZI_BUILD_TESTS)
  FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.9.0
  )
  FetchContent_MakeAvailable(catch2)
  include(CTest)
  include(Catch)

  add_executable(test_template tests/test_template.cpp)
  target_link_libraries(test_template PUBLIC Catch2::Catch2WithMain polyfit::polyfit)
  target_include_directories(test_template PRIVATE ${BAOBZI_INCLUDES})
  catch_discover_tests(test_template)

  if (BAOBZI_BUILD_SHARED)
    add_executable(test_c tests/test_c.cpp)
    target_link_libraries(test_c PUBLIC Catch2::Catch2WithMain baobzi)
    catch_discover_tests(test_c)
  endif()
endif()

if (${BAOBZI_BUILD_MATLAB})
  find_package(Matlab REQUIRED)
  matlab_add_mex(NAME baobzi_mex SRC src/matlab/baobzi_mex.cpp LINK_TO baobzi)
  target_include_directories(baobzi_mex PUBLIC ${BAOBZI_INCLUDES})
  install(FILES ${PROJECT_SOURCE_DIR}/src/matlab/baobzi.m
    DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/matlab)
  install(TARGETS baobzi_mex DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/matlab)
endif()


install(TARGETS ${INSTALL_TARGETS})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/licenses/Baobzi)

if (NOT DEFINED SKBUILD)
  install(FILES ${PROJECT_SOURCE_DIR}/src/python/__init__.py
    RENAME baobzi.py
    DESTINATION ${CMAKE_INSTALL_PREFIX}/share/baobzi/python)
endif()

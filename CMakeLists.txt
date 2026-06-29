# TODO: add some comment at the top
# TODO: change description and project name

# Preamble ####################################################################
#
cmake_minimum_required(VERSION 3.24.0)
project(fptools
    VERSION 0.1.0
    DESCRIPTION "something here"
    LANGUAGES CXX
)

# Load helper functions
include(${CMAKE_SOURCE_DIR}/cmake/FPToolsFunctions.cmake)

# Guard against in-source builds, which pollute the source tree and break
# clean rebuilds.
if (CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR 
        "In-source builds are not supported.\n"
        "Please create a build directory and run cmake from there:\n"
        "   cmake -B build\n"
        "   cmake --build build\n"
        "Remove CMakeCache.txt and CMakeFiles/ if they were created here."
    )
endif ()

# Set "Release" as the default build type rather than CMake's historic "Debug"
# default. Valid options are: Release, Debug, RelWithDebInfo, MinSizeRel.
fp_tools_set_default_build_type("Release")


# C++ Standard ################################################################
#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Output Directories ##########################################################
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


# Main Executable Target ######################################################
#
# Sources are listed explicitly rather than GLOB_RECURSE so CMake detects new 
# or removed files automatically on reconfigure.
#
add_executable(fptools
    src/fp_tools.cpp
    src/analysis/msd.cpp
    src/analysis/rdf.cpp
    src/analysis/mfpt.cpp
    src/analysis/corr.cpp
    src/io/config.cpp
    src/io/logger.cpp
    src/io/formats/lammps.cpp
    src/io/formats/xyz.cpp
    src/io/formats/binary.cpp
    src/io/formats/vasp.cpp
)

# TODO: fix this so it's not needed or do SOMETHING about it
add_compile_definitions(SOURCE_DIR="${CMAKE_SOURCE_DIR}")


# Headers: include/ for public headers, ext/ for vendored dependencies (e.g. TOML)
target_include_directories(fptools PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/analysis
    ${CMAKE_SOURCE_DIR}/include/core
    ${CMAKE_SOURCE_DIR}/include/io
    ${CMAKE_SOURCE_DIR}/include/io/formats
    ${CMAKE_SOURCE_DIR}/ext
)

# Link stdc++fs ONLY if using GCC version 9.0 or earlier
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(fptools PRIVATE stdc++fs)
endif()


# Dependencies ################################################################
#
# Each dependency is isolated in its own cmake/dependencies/ file. This keeps
# the root CMakeLists.txt readable and makes it easy to swap out or upgrade
# individual dependencies.
#
include(${CMAKE_SOURCE_DIR}/cmake/dependencies/FFTW.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/dependencies/OpenMP.cmake)

if (FPTOOLS_USE_MPI)
    include(${CMAKE_SOURCE_DIR}/cmake/dependencies/MPI.cmake)
endif ()


# Warnings ####################################################################
#
fp_tools_set_cxx_warnings(fptools)


# Tests #######################################################################
#
# GoogleTest is fetched automatically via FetchContent if not found on system.
# Tests are only built when FPTOOLS_BUILD_TESTS=ON (default).
# Run with: ctest --test-dir build
#
if (FPTOOLS_BUILD_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif ()


# Configuration Summary #######################################################
#
# Printed at the end of cmake configuration so build options are visible.
#
fp_tools_print_summary()

























# # CMake version requirement
# cmake_minimum_required(VERSION 3.16)
#
# # Project information
# project(fp-tools LANGUAGES CXX)
# set(FPTOOLS_VERSION_MAJOR 1)
# set(FPTOOLS_VERSION_MINOR 0)
#
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release CACHE STRING
#       "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
#       FORCE)
# endif()
#
# add_compile_definitions(SOURCE_DIR="${CMAKE_SOURCE_DIR}")
#
# # C++ standard
# set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)
#
# # ========= Configurable build options =========
# option(USE_MPI     "Enable MPI support"     OFF)
# option(USE_OMP     "Enable OpenMP support"  OFF)
# option(USE_CUDA    "Enable CUDA support"    OFF)
#
# # ========= Include directories =========
# include_directories(${CMAKE_SOURCE_DIR}/include)
# include_directories(${CMAKE_SOURCE_DIR}/ext)
#
# # ========= Source files and Executable Directory =========
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
# add_executable(fptools ${SOURCES})
#
# # ========= MPI =========
# if(USE_MPI)
#     find_package(MPI REQUIRED)
#     if(MPI_FOUND)
#         message(STATUS "MPI found. Enabling MPI support.")
#         add_definitions(-DUSE_MPI)
#         include_directories(${MPI_INCLUDE_PATH})
#         set(EXTRA_LIBS ${EXTRA_LIBS} MPI::MPI_CXX)
#     else()
#         message(FATAL_ERROR "USE_MPI is ON but MPI was not found.")
#     endif()
# endif()
#
# # ========= OpenMP =========
# # if(USE_OMP)
# #     find_package(OpenMP REQUIRED)
# #     if(OpenMP_FOUND)
# #         message(STATUS "OpenMP found. Enabling OpenMP support.")
# #         add_definitions(-DUSE_OMP)
# #         set(EXTRA_LIBS ${EXTRA_LIBS} OpenMP::OpenMP_CXX)
# #     else()
# #         message(FATAL_ERROR "USE_OMP is ON but OpenMP was not found.")
# #     endif()
# # endif()
#
# # if (USE_OMP)
# #   # 1. Locate the OpenMP package (needs libomp and a clang or gcc that supports -fopenmp)
# #   find_package(OpenMP REQUIRED)
# #
# #   # 2. Tell *each* target that uses OpenMP to link against it
# #   target_link_libraries(fptools PRIVATE OpenMP::OpenMP_CXX)
# #
# #   # 3. Pass the pre‑processor symbol so your #ifdef blocks are active
# #   target_compile_definitions(fptools PRIVATE USE_OMP)
# # endif()
#
# if(USE_OMP)
#   message(STATUS "Using OpenMP with libgomp (GNU OpenMP runtime)")
#
#   set(LIBGOMP_INCLUDE_DIR "/opt/homebrew/opt/gcc/include")
#   set(LIBGOMP_LIB_DIR "/opt/homebrew/Cellar/gcc/15.1.0/lib/gcc/15")
#
#   target_compile_options(fptools PRIVATE
#     -fopenmp=libgomp
#     -I${LIBGOMP_INCLUDE_DIR}
#   )
#
#   target_link_options(fptools PRIVATE
#     -fopenmp=libgomp
#     -L${LIBGOMP_LIB_DIR}
#     -lgomp
#   )
#
#   target_compile_definitions(fptools PRIVATE USE_OMP)
# endif()
#
#
#
# # ========= CUDA =========
# if(USE_CUDA)
#     find_package(CUDA REQUIRED)
#     if(CUDA_FOUND)
#         message(STATUS "CUDA found. Enabling CUDA support.")
#         add_definitions(-DUSE_CUDA)
#         include_directories(${CUDA_INCLUDE_DIRS})
#         set(EXTRA_LIBS ${EXTRA_LIBS} ${CUDA_LIBRARIES})
#     else()
#         message(FATAL_ERROR "USE_CUDA is ON but CUDA was not found.")
#     endif()
# endif()
#
# # ========= FFTW =========
#
# find_library(FFTW_LIB fftw3 REQUIRED)
# find_path(FFTW_INCLUDE_DIR fftw3.h REQUIRED)
#
# include_directories(${FFTW_INCLUDE_DIR})
# target_link_libraries(fptools ${FFTW_LIB})
#
# # find_package(FFTW REQUIRED COMPONENTS DOUBLE)
# # target_link_libraries(fptools PRIVATE FFTW::Double)
#
# # ========= Optimization flags (optional) =========
# target_compile_options(fptools PRIVATE -O3 -Wall -Wno-comment)

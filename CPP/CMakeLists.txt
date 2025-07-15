cmake_minimum_required(VERSION 3.18)
project(tic_tac LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL CXX)
find_package(BLAS REQUIRED)
find_library(OPENBLAS_LIB NAMES openblas REQUIRED)
find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED)
find_path(LAPACKE_INCLUDE_DIR lapacke.h PATHS /usr/include)

# Subdirectories
add_subdirectory(src/gsl_helpers)
add_subdirectory(src/blas_lapack_interface)
add_subdirectory(src)

# Executable
add_executable(tictac src/main.cpp)

# Link to the main static library
target_link_libraries(tictac PRIVATE tic_tac)

# Add include path for internal headers used directly in main.cpp
target_include_directories(tictac
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

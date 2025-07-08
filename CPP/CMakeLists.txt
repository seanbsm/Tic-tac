cmake_minimum_required(VERSION 3.18)
project(tic_tac LANGUAGES C CXX Fortran)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Collect source files
file(GLOB_RECURSE CPP_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE F90_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/*.f90)
file(GLOB_RECURSE F77_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/*.f)

# Create target
add_library(tic_tac STATIC ${CPP_SOURCES} ${F90_SOURCES} ${F77_SOURCES})

# Tell GSL to use external CBLAS
target_compile_definitions(tic_tac PRIVATE GSL_CBLAS_H) 

# Enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP")
    target_compile_options(tic_tac PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(tic_tac PRIVATE OpenMP::OpenMP_CXX)
endif()

# GSL (custom path)
set(GSL_DIR "$ENV{HOME}/local/gsl")
find_path(GSL_INCLUDE_DIR gsl/gsl_math.h PATHS ${GSL_DIR}/include)
find_library(GSL_LIBRARY NAMES gsl PATHS ${GSL_DIR}/lib)

message(STATUS "GSL include dir: ${GSL_INCLUDE_DIR}")
message(STATUS "GSL library: ${GSL_LIBRARY}")

# Optimized BLAS
find_package(BLAS REQUIRED)
message(STATUS "BLAS libraries found: ${BLAS_LIBRARIES}")

# LAPACK (LAPACKE)
find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED)
find_path(LAPACKE_INCLUDE_DIR lapacke.h PATHS /usr/include)
message(STATUS "LAPACKE library found: ${LAPACKE_LIBRARY}")
message(STATUS "LAPACKE include dir: ${LAPACKE_INCLUDE_DIR}")

# HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL CXX)
message(STATUS "HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")

# Include directories
target_include_directories(tic_tac
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${HDF5_INCLUDE_DIRS}
        ${LAPACKE_INCLUDE_DIR}
        ${GSL_INCLUDE_DIR}
)

# Compiler flags
target_compile_options(tic_tac PRIVATE
    -Wall -Wno-sign-compare -Wno-unused-variable -O3 -ggdb -fopenmp
)

# Fortran flags
set_source_files_properties(${F90_SOURCES} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC")
set_source_files_properties(${F77_SOURCES} PROPERTIES COMPILE_FLAGS "-O3")

# Link libraries
target_link_libraries(tic_tac
    PRIVATE
        ${BLAS_LIBRARIES}
        ${LAPACKE_LIBRARY}
        ${HDF5_LIBRARIES}
        ${GSL_LIBRARY}
        gfortran
        pthread
        m
        dl
)

# Optional install
install(TARGETS tic_tac DESTINATION lib)
install(DIRECTORY src/ DESTINATION include FILES_MATCHING PATTERN "*.h")

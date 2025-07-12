cmake_minimum_required(VERSION 3.18)
project(tic_tac LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------------
# GSL-only helper library
# -------------------------------
set(GSL_HELPER_SOURCES
    ${CMAKE_SOURCE_DIR}/General_functions/coupling_coefficients.cpp
    ${CMAKE_SOURCE_DIR}/auxiliary.cpp
)

add_library(gsl_helpers STATIC ${GSL_HELPER_SOURCES})

# GSL (no OpenBLAS)
set(GSL_DIR "$ENV{HOME}/local/gsl")
find_path(GSL_INCLUDE_DIR gsl/gsl_math.h PATHS ${GSL_DIR}/include)
find_library(GSL_LIBRARY NAMES gsl PATHS ${GSL_DIR}/lib)

target_include_directories(gsl_helpers
    PUBLIC ${GSL_INCLUDE_DIR}
)

target_link_libraries(gsl_helpers
    PRIVATE ${GSL_LIBRARY}
)

# -------------------------------
# Main tic_tac library
# -------------------------------
file(GLOB_RECURSE CPP_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/*.cpp)
file(GLOB_RECURSE F90_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/*.f90)
file(GLOB_RECURSE F77_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/*.f)

# Remove GSL-related files from main sources
list(REMOVE_ITEM CPP_SOURCES ${GSL_HELPER_SOURCES})

add_library(tic_tac STATIC ${CPP_SOURCES} ${F90_SOURCES} ${F77_SOURCES})

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_compile_options(tic_tac PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(tic_tac PRIVATE OpenMP::OpenMP_CXX)
endif()

# Optimized BLAS (OpenBLAS)
find_package(BLAS REQUIRED)
message(STATUS "BLAS libraries found: ${BLAS_LIBRARIES}")

# LAPACK
find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED)
find_path(LAPACKE_INCLUDE_DIR lapacke.h PATHS /usr/include)

# HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL CXX)

target_include_directories(tic_tac
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${HDF5_INCLUDE_DIRS}
        ${LAPACKE_INCLUDE_DIR}
)

target_compile_options(tic_tac PRIVATE
    -Wall -Wno-sign-compare -Wno-unused-variable -O3 -ggdb -fopenmp -Wl,--as-needed
)

set_source_files_properties(${F90_SOURCES} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC")
set_source_files_properties(${F77_SOURCES} PROPERTIES COMPILE_FLAGS "-O3")

target_link_libraries(tic_tac
    PRIVATE
        gsl_helpers
        ${BLAS_LIBRARIES}
        ${LAPACKE_LIBRARY}
        ${HDF5_LIBRARIES}
        gfortran
        pthread
        m
        dl
)

install(TARGETS tic_tac DESTINATION lib)
install(FILES
    ${CMAKE_SOURCE_DIR}/constants.h
    ${CMAKE_SOURCE_DIR}/dis_io_routines.h
    ${CMAKE_SOURCE_DIR}/error_management.h
    ${CMAKE_SOURCE_DIR}/make_permutation_matrix.h
    ${CMAKE_SOURCE_DIR}/make_potential_matrix.h
    ${CMAKE_SOURCE_DIR}/make_pw_symm_states.h
    ${CMAKE_SOURCE_DIR}/make_resolvent.h
    ${CMAKE_SOURCE_DIR}/make_swp_states.h
    ${CMAKE_SOURCE_DIR}/make_wp_states.h
    ${CMAKE_SOURCE_DIR}/run_organizer.h 
    ${CMAKE_SOURCE_DIR}/set_run_parameters.h 
    ${CMAKE_SOURCE_DIR}/set_run_parameters.h 
    ${CMAKE_SOURCE_DIR}/solve_faddeev.h 
    ${CMAKE_SOURCE_DIR}/type_defs.h 
    ${CMAKE_SOURCE_DIR}/General_functions/gauss_legendre.h 
    ${CMAKE_SOURCE_DIR}/General_functions/kinetic_conversion.h 
    ${CMAKE_SOURCE_DIR}/General_functions/legendre.h
    ${CMAKE_SOURCE_DIR}/General_functions/matrix_routines.h
    ${CMAKE_SOURCE_DIR}/General_functions/templates.h
    DESTINATION include
)


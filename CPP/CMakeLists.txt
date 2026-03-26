cmake_minimum_required(VERSION 3.18)
project(tic_tac LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found")

# Find HDF5 in a portable way
find_package(HDF5 REQUIRED COMPONENTS C HL CXX)

# Use imported targets if available, otherwise create a fallback target
if(TARGET HDF5::HDF5)
  message(STATUS "HDF5: using imported target HDF5::HDF5")
  set(HDF5_TARGET HDF5::HDF5)
elseif(TARGET hdf5::hdf5_cpp)
  message(STATUS "HDF5: using imported target hdf5::hdf5_cpp")
  set(HDF5_TARGET hdf5::hdf5_cpp)
else()
  message(STATUS "HDF5: using legacy variables fallback")
  add_library(HDF5::HDF5_CXX INTERFACE IMPORTED)
  target_include_directories(HDF5::HDF5_CXX INTERFACE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(HDF5::HDF5_CXX INTERFACE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  if(HDF5_LIBRARY_DIRS)
    target_link_directories(HDF5::HDF5_CXX INTERFACE ${HDF5_LIBRARY_DIRS})
  endif()
  set(HDF5_TARGET HDF5::HDF5_CXX)
endif()

# Check for HL target
if(TARGET hdf5::hdf5_hl)
  set(HDF5_HL_TARGET hdf5::hdf5_hl)
else()
  set(HDF5_HL_TARGET "")
endif()

find_package(BLAS REQUIRED)
message(STATUS "BLAS found")

find_library(OPENBLAS_LIB NAMES openblas REQUIRED)
if(OPENBLAS_LIB)
    message(STATUS "OpenBLAS found: ${OPENBLAS_LIB}")
else()
    message(STATUS "OpenBLAS not found")
endif()
find_library(LAPACKE_LIBRARY NAMES lapacke REQUIRED)
message(STATUS "LAPACKE found: ${LAPACKE_LIBRARY}")

find_path(LAPACKE_INCLUDE_DIR lapacke.h PATHS /usr/include)


# Add HDF5 includes globally (Debian/Ubuntu HDF5 1.10.x)
# Removed: now handled by target


# Create an interface target for OpenMP
add_library(OpenMPHelper INTERFACE)
target_link_libraries(OpenMPHelper INTERFACE OpenMP::OpenMP_CXX)
target_compile_options(OpenMPHelper INTERFACE ${OpenMP_CXX_FLAGS})


# Subdirectories
add_subdirectory(src/gsl_helpers)
add_subdirectory(src/blas_lapack_interface)
add_subdirectory(src)

# Executable
add_executable(tictac src/main.cpp)

# Link to the main static library
target_link_libraries(tictac PRIVATE tic_tac ${HDF5_TARGET} ${HDF5_HL_TARGET})

# Add include path for internal headers used directly in main.cpp
target_include_directories(tictac
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
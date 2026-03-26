# List all source files explicitly
set(UTILS_SOURCES
    disk_io_routines.cpp
    gauss_legendre.cpp
    kinetic_conversion.cpp
    legendre.cpp
    matrix_routines.cpp
)

# Create an object library from the collected sources
add_library(utils OBJECT ${UTILS_SOURCES})

# Set compile options
target_compile_options(utils PRIVATE -O3 -fPIC)

# Add private include path for internal headers (e.g., from other src/ subfolders)
target_include_directories(utils
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Link OpenMP flags and libraries via the interface target
if(TARGET HDF5::HDF5)
  target_link_libraries(utils PRIVATE OpenMPHelper HDF5::HDF5 ${HDF5_HL_TARGET})
elseif(TARGET hdf5::hdf5_cpp)
  target_link_libraries(utils PRIVATE OpenMPHelper hdf5::hdf5_cpp ${HDF5_HL_TARGET})
else()
  target_link_libraries(utils PRIVATE OpenMPHelper HDF5::HDF5_CXX)
endif()
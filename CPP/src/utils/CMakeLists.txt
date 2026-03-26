# Recursively collect all .cpp source files in this directory
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS *.cpp)

# Create an object library from the collected sources
add_library(utils OBJECT ${SOURCES})

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
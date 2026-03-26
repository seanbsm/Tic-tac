# Recursively collect all .cpp source files in this directory
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS *.cpp)

# Create an object library from the collected sources
add_library(core OBJECT ${SOURCES})

# Add private include path for internal headers (e.g., utils/constants.h)
target_include_directories(core
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Link OpenMP flags and libraries via the interface target
if(TARGET HDF5::HDF5)
  target_link_libraries(core PRIVATE OpenMPHelper HDF5::HDF5 ${HDF5_HL_TARGET})
elseif(TARGET hdf5::hdf5_cpp)
  target_link_libraries(core PRIVATE OpenMPHelper hdf5::hdf5_cpp ${HDF5_HL_TARGET})
else()
  target_link_libraries(core PRIVATE OpenMPHelper HDF5::HDF5_CXX)
endif()

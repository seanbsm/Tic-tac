# Recursively collect all .cpp source files in this directory
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS *.cpp)

# Create an object library from the collected sources
add_library(master OBJECT ${SOURCES})

# Add private include path for internal headers (e.g., from utils/, core/, etc.)
target_include_directories(master
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Link HDF5 for includes
if(TARGET HDF5::HDF5)
  target_link_libraries(master PRIVATE HDF5::HDF5 ${HDF5_HL_TARGET})
elseif(TARGET hdf5::hdf5_cpp)
  target_link_libraries(master PRIVATE hdf5::hdf5_cpp ${HDF5_HL_TARGET})
else()
  target_link_libraries(master PRIVATE HDF5::HDF5_CXX)
endif()

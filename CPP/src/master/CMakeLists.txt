# List all source files explicitly
set(MASTER_SOURCES
    error_management.cpp
    run_organizer.cpp
    set_run_parameters.cpp
)

# Create an object library from the collected sources
add_library(master OBJECT ${MASTER_SOURCES})

# Set compile options
target_compile_options(master PRIVATE -O3 -fPIC)

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

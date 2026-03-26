# List all source files explicitly
set(CORE_SOURCES
    make_permutation_matrix.cpp
    make_potential_matrix.cpp
    make_pw_symm_states.cpp
    make_resolvent.cpp
    make_swp_states.cpp
    make_wp_states.cpp
    solve_faddeev.cpp
)

# Create an object library from the collected sources
add_library(core OBJECT ${CORE_SOURCES})

# Set compile options
target_compile_options(core PRIVATE -O3 -fPIC)

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

# src/CMakeLists.txt

# Add subdirectories for each component
add_subdirectory(core)
add_subdirectory(master)
add_subdirectory(utils)
add_subdirectory(interactions)

# Create the main static library by linking object libraries
add_library(tic_tac STATIC)

# Link object libraries from subdirectories
target_link_libraries(tic_tac
    PRIVATE
        core
        master
        utils
        interactions
        gsl_helpers
        blas_lapack_interface
        ${HDF5_LIBRARIES}
        ${HDF5_HL_LIBRARIES}
        gfortran
        pthread
        m
        dl
)

# Include directories for HDF5 and LAPACKE
target_include_directories(tic_tac
    PRIVATE
        ${HDF5_INCLUDE_DIRS}
        ${LAPACKE_INCLUDE_DIR}
)

# Enable OpenMP if available
if(OpenMP_CXX_FOUND)
    target_compile_options(tic_tac PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(tic_tac PRIVATE OpenMP::OpenMP_CXX)
endif()

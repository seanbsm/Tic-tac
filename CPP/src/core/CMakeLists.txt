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
target_link_libraries(core PRIVATE OpenMPHelper)

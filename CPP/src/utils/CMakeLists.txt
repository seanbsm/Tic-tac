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
target_link_libraries(utils PRIVATE OpenMPHelper)
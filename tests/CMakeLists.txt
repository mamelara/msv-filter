# tests/CMakeLists.txt
# Google Test setup using FetchContent (no find_package)

include(FetchContent)

# Fetch Google Test from GitHub
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2  # Use a specific version for reproducibility
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Make Google Test available
FetchContent_MakeAvailable(googletest)

# Create test executable
add_executable(msv_tests
    test_msv_basic.cpp
    test_msv_edge_cases.cpp
    stub_msv.cpp  # Stub implementation - user will replace with actual algorithm
)

# Link against Google Test
target_link_libraries(msv_tests
    GTest::gtest
    GTest::gtest_main
)

# Include directories
target_include_directories(msv_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# Add additional source files from main project that tests depend on
target_sources(msv_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src/aa_alphabet.cpp
)

# Discover and register tests with CTest
include(GoogleTest)
gtest_discover_tests(msv_tests)

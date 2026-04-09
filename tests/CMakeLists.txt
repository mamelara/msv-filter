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

if(TARGET msv_hmmer_bridge)
    add_executable(msv_tests
        test_msv_hmmer_integration.cpp
        ${CMAKE_SOURCE_DIR}/src/aa_alphabet.cpp
        ${CMAKE_SOURCE_DIR}/src/msv_filter.cpp
    )

    target_link_libraries(msv_tests
        GTest::gtest
        GTest::gtest_main
        msv_hmmer_bridge
    )

    target_include_directories(msv_tests PRIVATE
        ${CMAKE_SOURCE_DIR}/include
    )

    include(GoogleTest)
    gtest_discover_tests(msv_tests)
else()
    message(WARNING "Skipping tests: msv_hmmer_bridge target is required for integration tests")
endif()

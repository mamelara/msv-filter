cmake_minimum_required(VERSION 4.0)
project(msv_filter)

set(CMAKE_CXX_STANDARD 17)

# Enable clang-tidy if available
find_program(CLANG_TIDY_EXE NAMES clang-tidy)
if(CLANG_TIDY_EXE)
    message(STATUS "Found clang-tidy: ${CLANG_TIDY_EXE}")
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
else()
    message(STATUS "clang-tidy not found")
endif()

# Main executable
add_executable(msv_filter
        src/main.cpp
        src/aa_alphabet.cpp
)

target_include_directories(msv_filter PRIVATE include)

# Enable testing and add tests subdirectory
enable_testing()
add_subdirectory(tests)

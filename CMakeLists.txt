cmake_minimum_required(VERSION 4.0)
project(msv_filter)

set(CMAKE_CXX_STANDARD 17)

# Main executable
add_executable(msv_filter
        src/main.cpp
        src/aa_alphabet.cpp
)

target_include_directories(msv_filter PRIVATE include)

# Enable testing and add tests subdirectory
enable_testing()
add_subdirectory(tests)

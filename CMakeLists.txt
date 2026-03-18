cmake_minimum_required(VERSION 4.0)
project(msv_filter LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)

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
        src/msv_filter.cpp
        src/aa_alphabet.cpp
)

target_include_directories(msv_filter PRIVATE include)

option(MSV_ENABLE_HMMER_BRIDGE "Enable HMMER/Easel bridge integration" OFF)

if(MSV_ENABLE_HMMER_BRIDGE)
    set(MSV_HMMER_PREFIX "" CACHE PATH "Install prefix for HMMER/Easel (must contain include/ and lib/)")

    find_path(MSV_HMMER_INCLUDE_DIR
        NAMES hmmer.h
        HINTS
            ${MSV_HMMER_PREFIX}
        PATH_SUFFIXES include
    )

    find_library(MSV_HMMER_LIBRARY
        NAMES hmmer
        HINTS
            ${MSV_HMMER_PREFIX}
        PATH_SUFFIXES lib
    )

    find_library(MSV_EASEL_LIBRARY
        NAMES easel
        HINTS
            ${MSV_HMMER_PREFIX}
        PATH_SUFFIXES lib
    )

    if(NOT MSV_HMMER_INCLUDE_DIR OR NOT MSV_HMMER_LIBRARY OR NOT MSV_EASEL_LIBRARY)
        message(FATAL_ERROR
            "MSV_ENABLE_HMMER_BRIDGE is ON, but HMMER/Easel were not found. "
            "Set MSV_HMMER_PREFIX to a valid install prefix (for example from ./configure --prefix=...).")
    endif()

    add_library(msv_hmmer_bridge STATIC
        src/hmmer_c_bridge.c
        src/msv_bridge_adapter.cpp
    )

    target_include_directories(msv_hmmer_bridge
        PUBLIC
            include
        PRIVATE
            ${MSV_HMMER_INCLUDE_DIR}
    )

    target_link_libraries(msv_hmmer_bridge
        PUBLIC
            ${MSV_HMMER_LIBRARY}
            ${MSV_EASEL_LIBRARY}
    )

    add_executable(hmmer_parse_upstream
        src/hmmer_parse_upstream.cpp
        src/aa_alphabet.cpp
    )

    target_include_directories(hmmer_parse_upstream
        PRIVATE
            include
    )

    target_link_libraries(hmmer_parse_upstream
        PRIVATE
            msv_hmmer_bridge
    )
endif()

option(MSV_ENABLE_TESTS "Enable unit tests" ON)

if(MSV_ENABLE_TESTS)
    enable_testing()
    if(EXISTS "${CMAKE_SOURCE_DIR}/tests/CMakeLists.txt")
        add_subdirectory(tests)
    else()
        message(WARNING "Skipping tests: tests/CMakeLists.txt was not found in working tree")
    endif()
endif()

cmake_minimum_required(VERSION 4.0)

if(NOT DEFINED CMAKE_C_COMPILER)
    find_program(MSV_GCC_15 NAMES gcc-15)
    if(MSV_GCC_15)
        set(CMAKE_C_COMPILER "${MSV_GCC_15}" CACHE FILEPATH "Preferred C compiler" FORCE)
    endif()
endif()

if(NOT DEFINED CMAKE_CXX_COMPILER)
    find_program(MSV_GXX_15 NAMES g++-15)
    if(MSV_GXX_15)
        set(CMAKE_CXX_COMPILER "${MSV_GXX_15}" CACHE FILEPATH "Preferred CXX compiler" FORCE)
    endif()
endif()

project(msv_filter LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)

find_package(OpenMP REQUIRED)

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
    set(MSV_HMMER_PREFIX "" CACHE PATH "HMMER root: install prefix (include/lib) or source tree root (src/easel)")

    find_path(MSV_HMMER_INCLUDE_DIR
        NAMES hmmer.h
        HINTS
            ${MSV_HMMER_PREFIX}
        PATH_SUFFIXES include src
    )

    find_path(MSV_EASEL_INCLUDE_DIR
        NAMES easel.h
        HINTS
            ${MSV_HMMER_PREFIX}
        PATH_SUFFIXES include easel
    )

    find_library(MSV_HMMER_LIBRARY
        NAMES hmmer libhmmer.a
        HINTS
            ${MSV_HMMER_PREFIX}
        PATH_SUFFIXES lib src
    )

    find_library(MSV_EASEL_LIBRARY
        NAMES easel libeasel.a
        HINTS
            ${MSV_HMMER_PREFIX}
        PATH_SUFFIXES lib easel
    )

    if(NOT MSV_HMMER_INCLUDE_DIR OR NOT MSV_EASEL_INCLUDE_DIR OR NOT MSV_HMMER_LIBRARY OR NOT MSV_EASEL_LIBRARY)
        message(FATAL_ERROR
            "MSV_ENABLE_HMMER_BRIDGE is ON, but HMMER/Easel were not found. "
            "Set MSV_HMMER_PREFIX to either: "
            "(a) an install prefix from ./configure --prefix=..., or "
            "(b) a built source tree root containing src/libhmmer.a and easel/libeasel.a.")
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
            ${MSV_EASEL_INCLUDE_DIR}
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

    add_executable(msv_filter_scan
        src/msv_filter_scan.cpp
        src/aa_alphabet.cpp
        src/msv_filter.cpp
    )

    target_include_directories(msv_filter_scan
        PRIVATE
            include
    )

    target_link_libraries(msv_filter_scan
        PRIVATE
            msv_hmmer_bridge
            OpenMP::OpenMP_CXX
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

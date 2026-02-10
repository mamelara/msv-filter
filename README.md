# msv-filter

A C++ implementation of the MSV (Multiple Segment Viterbi) filter algorithm for HMMER (Hidden Markov Model for protein sequence analysis).

## Overview

This project provides a re-implementation of the `p7_GMSV` function from HMMER, which performs MSV (Multiple Segment Viterbi) scoring for biological sequence alignment. The MSV filter is a computationally efficient first-pass filter used in HMMER's protein homology search pipeline.

### Key Components

- **HMMER-compatible types** (`hmmer_types.hpp`): Replicates essential structures from HMMER
- **Amino acid alphabet** (`aa_alphabet.cpp/hpp`): Digital sequence encoding
- **Profile handling** (`profile.hpp`): HMM profile structure management
- **DP Matrix** (`dp_matrix.hpp`): Dynamic programming matrix for alignment scoring
- **Mock data generator** (`mock_data.hpp`): Test data generation utilities

## Building the Project

### Prerequisites

- CMake 4.0 or higher
- C++17 compatible compiler
- Internet connection (for fetching Google Test)

### Build Instructions

1. **Create a build directory:**
```bash
mkdir cmake-build-test
cd cmake-build-test
```

2. **Configure with CMake:**
```bash
cmake ..
```

3. **Build the project:**
```bash
cmake --build .
```

This will build:
- `msv_filter` - The main executable demonstrating mock inputs
- `msv_tests` - The unit test executable (uses Google Test)

### Build Output

After building, you'll find:
- Main executable: `cmake-build-test/msv_filter`
- Test executable: `cmake-build-test/tests/msv_tests`

## Running the Application

### Run the Main Program

```bash
./cmake-build-test/msv_filter
```

This demonstrates the creation of mock inputs for the `p7_GMSV` function, including:
- Digital sequence generation
- HMM profile creation
- DP matrix allocation
- Memory layout visualization

## Running Tests

### Using CTest (Recommended)

From the build directory:
```bash
cd cmake-build-test
ctest --output-on-failure
```

### Running Tests Directly

```bash
./cmake-build-test/tests/msv_tests
```

### Test Output

The test suite includes 27 tests covering:
- **Basic functionality**: Constant patterns, single position models, alternating patterns
- **Edge cases**: Empty sequences, large scores, boundary conditions
- **Verification tests**: Test vector validation, sentinel preservation

Example output:
```
Test project /path/to/msv-filter/cmake-build-test
      Start  1: MSVBasicTest.ConstantAllOnes
 1/27 Test  #1: MSVBasicTest.ConstantAllOnes ..........................   Passed    0.01 sec
...
89% tests passed, 3 tests failed out of 27
```

Note: Some tests may fail if the MSV algorithm implementation is not yet complete (currently using stub implementation).

## CLion Setup

### Opening the Project

1. Open CLion
2. Select **File → Open**
3. Navigate to the project root directory (`msv-filter`)
4. Click **OK**

### CMake Configuration

CLion will automatically detect the `CMakeLists.txt` and configure the project. The default build directory will be `cmake-build-debug` or `cmake-build-release`.

### Running in CLion

**Main Application:**
1. In the top toolbar, select `msv_filter` from the run configuration dropdown
2. Click the **Run** button (green triangle) or press `Ctrl+R` (macOS) / `Shift+F10` (Linux/Windows)

**Tests:**
1. Select `msv_tests` from the run configuration dropdown
2. Click **Run** to execute all tests
3. Or use **Tools → CMake → Run Tests** to run via CTest

### Debugging

1. Set breakpoints in the code by clicking in the gutter
2. Select the target (`msv_filter` or `msv_tests`) from the run configuration
3. Click the **Debug** button (bug icon) or press `Ctrl+D` (macOS) / `Shift+F9` (Linux/Windows)

### CMake Options in CLion

To customize CMake options:
1. Go to **File → Settings** (or **CLion → Preferences** on macOS)
2. Navigate to **Build, Execution, Deployment → CMake**
3. Add options in **CMake options** field if needed

## Code Quality Tools

This project uses `clang-tidy` and `clang-format` for code quality and consistency.

### Prerequisites

Install LLVM/Clang tools:
- **macOS**: `brew install llvm` (includes clang-tidy and clang-format)
- **Ubuntu/Debian**: `sudo apt-get install clang-tidy clang-format`
- **Fedora**: `sudo dnf install clang-tools-extra`

### clang-tidy

**Automatic Integration:**
CMake automatically detects and enables clang-tidy during the build process. Warnings will appear during compilation.

**Manual Usage:**
```bash
# Run on all source files
clang-tidy src/*.cpp include/*.hpp tests/*.cpp -- -Iinclude

# Run with automatic fixes (use with caution)
clang-tidy src/*.cpp --fix -- -Iinclude

# Run specific checks only
clang-tidy src/*.cpp -checks='cppcoreguidelines-*,modernize-*' -- -Iinclude
```

**Configuration:**
- Config file: `.clang-tidy`
- Enables most checks by default
- Excludes vendor-specific and overly strict checks
- Configured for C++17

### clang-format

**Format all source files:**
```bash
# Format in place
clang-format -i src/*.cpp include/*.hpp tests/*.cpp

# Check if files are formatted (CI-friendly)
clang-format --dry-run --Werror src/*.cpp include/*.hpp tests/*.cpp
```

**Configuration:**
- Config file: `.clang-format`
- Based on LLVM style
- 4-space indentation
- 120 character line limit
- Right-aligned pointers

### CLion Integration

Both tools are automatically integrated with CLion:

**clang-tidy:**
- Warnings appear as you type
- Navigate to **Settings → Editor → Inspections → C/C++ → Clang-Tidy**

**clang-format:**
- Use **Code → Reformat Code** (Cmd+Option+L on macOS, Ctrl+Alt+L on Linux/Windows)
- Configure in **Settings → Editor → Code Style → C/C++ → Gear Icon → Import Scheme → .clang-format file**

## Project Structure

```
msv-filter/
├── CMakeLists.txt          # Main CMake configuration
├── README.md               # This file
├── src/                    # Source files
│   ├── main.cpp           # Main executable
│   └── aa_alphabet.cpp    # Amino acid alphabet implementation
├── include/               # Header files
│   ├── hmmer_types.hpp    # HMMER-compatible type definitions
│   ├── aa_alphabet.hpp    # Alphabet definitions
│   ├── profile.hpp        # Profile structures
│   ├── dp_matrix.hpp      # DP matrix implementation
│   └── mock_data.hpp      # Mock data generation
└── tests/                 # Unit tests
    ├── CMakeLists.txt     # Test CMake configuration
    ├── test_msv_basic.cpp # Basic functionality tests
    ├── test_msv_edge_cases.cpp # Edge case tests
    └── stub_msv.cpp       # Stub MSV implementation
```

## HMMER Integration

This project is designed to be compatible with HMMER's data structures:

- `DigitalResidue` - Digital sequence encoding (uint8_t)
- `P7_PROFILE` - Profile structure with match scores and transitions
- `P7_GMX` - Generic DP matrix for alignment
- MSV algorithm - Uses match scores only (MSC), ignores transitions

## License

This project is designed for educational and research purposes related to HMMER protein sequence analysis.

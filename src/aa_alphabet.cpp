#include "aa_alphabet.hpp"

AminoAcidAlphabet::AminoAcidAlphabet() {
    // --- A. Define Inputs ---
    // As passed in the function call: (..., "ACDEF...", 20, 29)
    const char* alphabet_str = "ACDEFGHIKLMNPQRSTVWY-BJZOUX*~";
    K = 20;
    Kp = 29;
    type = eslNONSTANDARD;
    sym = alphabet_str;

    // --- B. Allocation ---
    inmap.resize(128, digitalResidueIllegal); // Default everything to ILLEGAL
    ndegen.resize(Kp, 0);              // Default count to 0
    degen.resize(Kp * K, 0);           // Default matrix to 0

    // --- C. Initialize Input Map (ASCII -> Index) ---
    // Only maps the exact characters in the string (Case sensitive per the C code)
    for (int x = 0; x < Kp; x++) {
      char c = sym[x];
      inmap[(int)c] = x;
    }

    // --- D. Initialize Degeneracy Logic ---

    // 1. Base Alphabet (0..19)
    // Maps uniquely to itself.
    for (int x = 0; x < K; x++) {
      ndegen[x] = 1;
      set_degen(x, x, 1); // degen[x][x] = 1
    }

    // 2. The "Any" Character Logic
    // The code defines 'any' as index (Kp - 3).
    // For size 29, Kp-3 = 26.
    // In string "ACDEFGHIKLMNPQRSTVWY-BJZOUX*~", index 26 is 'X'.
    int any_idx = Kp - 3;

    ndegen[any_idx] = K; // X represents all 20 amino acids
    for (int y = 0; y < K; y++) {
      set_degen(any_idx, y, 1); // X matches A, C, D, E...
    }

    // Note: The C function leaves B, J, Z, etc. with ndegen=0.
    // It does not apply biological logic (like B=D|N) because it's the "Custom" function.
};

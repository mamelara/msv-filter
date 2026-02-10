#ifndef MSV_FILTER_AA_ALPHABET_H
#define MSV_FILTER_AA_ALPHABET_H
#include <string>
#include <vector>
#include "hmmer_types.hpp"
const int eslNONSTANDARD = 0;   // Type tag (generic)

struct AminoAcidAlphabet {
  // 1. Core Dimensions
  int K;  // Canonical size (20)
  int Kp; // Total size (29)
  int type;

  // 2. Data Structures
  std::string sym;             // The symbol string "ACDEF..."
  std::vector<int> inmap;      // Maps ASCII (0-127) -> Digital Index
  std::vector<int> ndegen;     // How many residues does this char represent?
  std::vector<uint8_t> degen;  // The degeneracy matrix (Flattened 2D: Kp rows * K cols)

  AminoAcidAlphabet();

  void set_degen(int row, int col, uint8_t val) {
    degen[(row * K) + col] = val;
  }

  uint8_t get_degen(int row, int col) const {
    return degen[(row * K) + col];
  }
};

#endif //MSV_FILTER_AA_ALPHABET_H

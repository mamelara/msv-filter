/*******************************************************************************
 * File: src/main.cpp
 * Description: Demonstrates mocked inputs for p7_GMSV function.
 * 
 * This file creates all the necessary inputs to call a re-implemented
 * p7_GMSV function with mocked data:
 *   - DigitalResidue digital sequence
 *   - P7_PROFILE structure
 *   - P7_GMX DP matrix
 *   - expected_hit_count parameter
 *   - msv_score return value
 ******************************************************************************/

#include <iostream>
#include <vector>
#include <cmath>
#include "hmmer_types.hpp"
#include "aa_alphabet.hpp"
#include "profile.hpp"
#include "dp_matrix.hpp"
#include "mock_data.hpp"

/*******************************************************************************
 * Example signature of the MSV function to be implemented:
 * 
 * int p7_GMSV(const DigitalResidue *digital_sequence, int sequence_length, const P7_PROFILE *gm, 
 *             P7_GMX *gx, float expected_hit_count, float *msv_score)
 * 
 * Where:
 *   - digital_sequence: Digital sequence (1-indexed array of residue indices)
 *   - sequence_length:   Sequence length
 *   - gm:  Profile with match scores (MSV only uses MSC, not transitions)
 *   - gx:  DP matrix with match states and special states
 *   - expected_hit_count:  Expected number of hits (typically 2.0)
 *   - msv_score: Return value for MSV score
 ******************************************************************************/

float compute_msv(const DigitalResidue* digital_sequence,
  int sequence_length, const HMMProfile& profile, DPMatrix& dp_matrix, const float expected_hit_count, const float msv_score) {

  float loop_transition_score = logf(static_cast<float>(sequence_length) / static_cast<float>(sequence_length + 3));
  float move_transition_score = logf( 3.0f / static_cast<float>(sequence_length + 3));
  float begin_transition_score = logf( 2.0f / static_cast<float>(profile.model_length * profile.model_length + 1));
  float t_state_to_e_state = logf((expected_hit_count - 1.0f) / expected_hit_count);
  float e_state_to_c_state = logf(1.0f / expected_hit_count);

  dp_matrix.special(0, p7G_N) = 0.0f;
  dp_matrix.special(0, p7G_B) = move_transition_score;
  dp_matrix.special(0, p7G_E) = dp_matrix.special(0, p7G_C) = dp_matrix.special(0, p7G_J) = -eslINFINITY;
  for (int i = 0; i <= profile.model_length; i++) dp_matrix.special(0, i) = -eslINFINITY;

  for (int i = 1; i <= sequence_length; i++) {
    dp_matrix.match(i, 0) = -eslINFINITY;
    dp_matrix.special(i, p7G_E) = -eslINFINITY;
    for (int k = 1; k <= profile.model_length; k++) {
      DigitalResidue residue = digital_sequence[i];
      float msc = profile.match_score(k, residue);
      
      // MMX(i-1,k-1): previous match state
      float prev_match = dp_matrix.match(i-1, k-1);
      
      // XMX(i-1,p7G_B) + tbmk: from Begin state
      float from_begin = dp_matrix.special(i-1, p7G_B) + begin_transition_score;
      
      // ESL_MAX: take the maximum
      float best_prev = std::max(prev_match, from_begin);
      
      // MMX(i,k) = MSC(k) + max(...)
      dp_matrix.match(i, k) = msc + best_prev;
    }

    dp_matrix.special(i, p7G_J) = std::max((dp_matrix.special(i - 1, p7G_J) + loop_transition_score) ,
      (dp_matrix.special(i - 1, p7G_E) + t_state_to_e_state));

    dp_matrix.special(i, p7G_C) = std::max(dp_matrix.special(i, p7G_C) ,
      (dp_matrix.special(i, p7G_E) + e_state_to_c_state));

    dp_matrix.special(i, p7G_N) = dp_matrix.special(i - 1, p7G_N) + loop_transition_score;

    dp_matrix.special(i, p7G_B) = std::max((dp_matrix.special(i, p7G_N) + move_transition_score),
      (dp_matrix.special(i, p7G_J) + move_transition_score));
  }

}

int main() {

    AminoAcidAlphabet abc;

    int sequence_length = 15;  // Sequence length
    std::vector<DigitalResidue> digital_sequence = MockDataGenerator::create_simple_sequence(sequence_length, abc);
    MockDataGenerator::print_sequence(digital_sequence, sequence_length, abc);

    // --- Step 3: Create Mock Profile ---
    int model_length = 10;  // Model length
    HMMProfile profile = MockDataGenerator::create_simple_profile(model_length, abc);

    DPMatrix dp_matrix(model_length, sequence_length);

    float expected_hit_count = 2.0f;  // Expected number of hits
    float msv_score = 0.0f;  // Will hold the result

  msv_filter(digital_sequence.data(), sequence_length, profile, dp_matrix, expected_hit_count, msv_score);



    return 0;
}

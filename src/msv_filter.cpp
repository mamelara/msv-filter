#include <algorithm>
#include <cmath>
#include <iostream>

#include "msv_filter.hpp"

int msv_filter(const DigitalResidue* digital_sequence,
  int sequence_length,
  const HMMProfile& profile,
  DPMatrix& dp_matrix,
  const float expected_hit_count,
  float* out_msv_score) {
  if (out_msv_score == nullptr) {
    return -1;
  }

  float loop_transition_score = logf(static_cast<float>(sequence_length) / static_cast<float>(sequence_length + 3));
  float move_transition_score = logf( 3.0f / static_cast<float>(sequence_length + 3));
  float begin_transition_score = logf( 2.0f / static_cast<float>(profile.model_length * (profile.model_length + 1)));
  float e_state_to_j_state = logf((expected_hit_count - 1.0f) / expected_hit_count);
  float e_state_to_c_state = logf(1.0f / expected_hit_count);

  dp_matrix.special(0, p7G_N) = 0.0f;
  dp_matrix.special(0, p7G_B) = move_transition_score;
  dp_matrix.special(0, p7G_E) = -eslINFINITY;
  dp_matrix.special(0, p7G_C) = -eslINFINITY;
  dp_matrix.special(0, p7G_J) = -eslINFINITY;

  for (int i = 0; i <= profile.model_length; i++) dp_matrix.match(0, i) = -eslINFINITY;

  for (int i = 1; i <= sequence_length; i++) {
    DigitalResidue residue = digital_sequence[i];
    dp_matrix.match(i, 0) = -eslINFINITY;
    dp_matrix.special(i, p7G_E) = -eslINFINITY;
    for (int k = 1; k <= profile.model_length; k++) {
      float msc = profile.match_score(k, residue);

      // MMX(i-1,k-1): previous match state
      float prev_match = dp_matrix.match(i-1, k-1);

      // XMX(i-1,p7G_B) + tbmk: from Begin state
      float from_begin = dp_matrix.special(i-1, p7G_B) + begin_transition_score;

      // ESL_MAX: take the maximum
      float best_prev = std::max(prev_match, from_begin);

      // MMX(i,k) = MSC(k) + max(...)
      dp_matrix.match(i, k) = msc + best_prev;
      dp_matrix.special(i, p7G_E) = std::max(dp_matrix.special(i, p7G_E), dp_matrix.match(i, k));
    }

    dp_matrix.special(i, p7G_J) = std::max(dp_matrix.special(i-1, p7G_J) + loop_transition_score,
      dp_matrix.special(i, p7G_E) + e_state_to_j_state);

    dp_matrix.special(i, p7G_C) = std::max(dp_matrix.special(i-1, p7G_C) + loop_transition_score,
      (dp_matrix.special(i, p7G_E) + e_state_to_c_state));

    dp_matrix.special(i, p7G_N) = dp_matrix.special(i-1, p7G_N) + loop_transition_score;

    dp_matrix.special(i, p7G_B) = std::max((dp_matrix.special(i, p7G_N) + move_transition_score),
      (dp_matrix.special(i, p7G_J) + move_transition_score));
  }

  *out_msv_score = dp_matrix.special(sequence_length, p7G_C) + move_transition_score;

  return 0;
}

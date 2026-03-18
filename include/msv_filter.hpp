#ifndef MSV_FILTER_MSV_FILTER_HPP
#define MSV_FILTER_MSV_FILTER_HPP

#include "dp_matrix.hpp"
#include "hmmer_types.hpp"
#include "profile.hpp"

int msv_filter(const DigitalResidue* digital_sequence,
               int sequence_length,
               const HMMProfile& profile,
               DPMatrix& dp_matrix,
               float expected_hit_count,
               float* out_msv_score);

#endif

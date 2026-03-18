#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "aa_alphabet.hpp"
#include "hmmer_c_bridge.h"
#include "msv_bridge_adapter.hpp"
#include "msv_filter.hpp"

namespace {

const char* get_env_or_default(const char* name, const char* default_value) {
  const char* value = std::getenv(name);
  return (value == nullptr || value[0] == '\0') ? default_value : value;
}

int get_env_int_or_default(const char* name, int default_value) {
  const char* value = std::getenv(name);
  if (value == nullptr || value[0] == '\0') {
    return default_value;
  }

  return std::atoi(value);
}

float get_env_float_or_default(const char* name, float default_value) {
  const char* value = std::getenv(name);
  if (value == nullptr || value[0] == '\0') {
    return default_value;
  }

  return std::strtof(value, nullptr);
}

bool both_same_infinity(float lhs, float rhs) {
  return std::isinf(lhs) && std::isinf(rhs) && (std::signbit(lhs) == std::signbit(rhs));
}

std::string label_for_case(const char* hmm_name, const char* seq_name, int hmm_idx, int seq_idx) {
  return "hmm#" + std::to_string(hmm_idx) + "(" + hmm_name + ")" +
         " vs seq#" + std::to_string(seq_idx) + "(" + seq_name + ")";
}

}  // namespace

TEST(MSVHmmerIntegrationTest, QuickParityAgainstHmmerGMSV) {
  const char* hmm_path = get_env_or_default("MSV_ITEST_HMM_PATH", "databases/panther.100.hmm");
  const char* fasta_path = get_env_or_default("MSV_ITEST_FASTA_PATH", "inputs/Arabidopsis_thaliana.100.pep.fa");
  const int max_hmms = get_env_int_or_default("MSV_ITEST_MAX_HMMS", 3);
  const int max_seqs = get_env_int_or_default("MSV_ITEST_MAX_SEQS", 25);
  const float expected_hit_count = get_env_float_or_default("MSV_ITEST_NU", 2.0f);
  const float score_tolerance = get_env_float_or_default("MSV_ITEST_TOL", 1e-4f);

  ASSERT_GT(max_hmms, 0) << "MSV_ITEST_MAX_HMMS must be > 0";
  ASSERT_GT(max_seqs, 0) << "MSV_ITEST_MAX_SEQS must be > 0";

  MsvHmmReader* hmm_reader = nullptr;
  const int open_status = msv_hmm_reader_open(hmm_path, &hmm_reader);
  if (open_status != MSV_BRIDGE_OK) {
    GTEST_SKIP() << "Skipping integration test. Unable to open HMM path: " << hmm_path;
  }

  int processed_hmms = 0;
  int processed_pairs = 0;
  int within_tolerance_pairs = 0;
  int nonfinite_equal_pairs = 0;
  int mismatched_pairs = 0;
  float max_abs_error = 0.0f;
  float mean_abs_error = 0.0f;
  float error_sum = 0.0f;
  int finite_compare_pairs = 0;
  int finite_mismatch_pairs = 0;
  int nonfinite_mismatch_pairs = 0;
  std::string worst_case_label;
  std::vector<std::string> mismatch_samples;

  for (int hmm_idx = 0; hmm_idx < max_hmms; ++hmm_idx) {
    MsvHmmRecord* hmm_record = nullptr;
    const int next_hmm_status = msv_hmm_reader_next(hmm_reader, &hmm_record);
    if (next_hmm_status == MSV_BRIDGE_EOF) {
      break;
    }

    ASSERT_EQ(MSV_BRIDGE_OK, next_hmm_status)
        << "Failed reading HMM record: " << msv_hmm_reader_last_error(hmm_reader);
    ASSERT_NE(nullptr, hmm_record);
    ++processed_hmms;

    MsvSeqReader* seq_reader = nullptr;
    const int seq_open_status = msv_seq_reader_open(fasta_path, hmm_record, &seq_reader);
    ASSERT_EQ(MSV_BRIDGE_OK, seq_open_status)
        << "Failed opening FASTA path: " << fasta_path;

    for (int seq_idx = 0; seq_idx < max_seqs; ++seq_idx) {
      MsvSeqRecord* seq_record = nullptr;
      const int next_seq_status = msv_seq_reader_next(seq_reader, &seq_record);
      if (next_seq_status == MSV_BRIDGE_EOF) {
        break;
      }

      ASSERT_EQ(MSV_BRIDGE_OK, next_seq_status)
          << "Failed reading sequence: " << msv_seq_reader_last_error(seq_reader);
      ASSERT_NE(nullptr, seq_record);

      const std::string case_label = label_for_case(msv_hmmrecord_name(hmm_record),
                                                    msv_seqrecord_name(seq_record),
                                                    hmm_idx,
                                                    seq_idx);

      MsvProfileCtx* profile_ctx = nullptr;
      const int seq_length = static_cast<int>(msv_seqrecord_length(seq_record));
      const int ctx_status = msv_profilectx_create(hmm_record, seq_length, MSV_PROFILE_MODE_LOCAL, &profile_ctx);
      ASSERT_EQ(MSV_BRIDGE_OK, ctx_status) << "Failed creating profile context for " << case_label;
      ASSERT_NE(nullptr, profile_ctx);

      float hmmer_score = 0.0f;
      const int hmmer_status = msv_profilectx_score_gmsv(profile_ctx, seq_record, expected_hit_count, &hmmer_score);
      ASSERT_EQ(MSV_BRIDGE_OK, hmmer_status)
          << "generic p7_GMSV (hmmer/src/generic_msv.c) failed for " << case_label << ": "
          << msv_profilectx_last_error(profile_ctx);

      AminoAcidAlphabet alphabet;
      const std::vector<DigitalResidue> digital_sequence = make_digital_sequence_from_bridge(seq_record);
      HMMProfile cpp_profile = make_profile_from_bridge(profile_ctx, alphabet);
      DPMatrix cpp_dp = make_dp_matrix_from_bridge(seq_record, profile_ctx);

      float our_score = 0.0f;
      const int our_status = msv_filter(
          digital_sequence.data(),
          seq_length,
          cpp_profile,
          cpp_dp,
          expected_hit_count,
          &our_score);

      if (our_status != 0) {
        ++mismatched_pairs;
        if (mismatch_samples.size() < 5) {
          mismatch_samples.push_back(case_label + " status mismatch: ours=" + std::to_string(our_status));
        }
      }

      if (std::isfinite(hmmer_score) && std::isfinite(our_score)) {
        const float abs_error = std::fabs(hmmer_score - our_score);
        error_sum += abs_error;
        ++finite_compare_pairs;
        if (abs_error > max_abs_error) {
          max_abs_error = abs_error;
          worst_case_label = case_label;
        }
        if (abs_error <= score_tolerance) {
          ++within_tolerance_pairs;
        } else {
          ++mismatched_pairs;
          ++finite_mismatch_pairs;
          if (mismatch_samples.size() < 5) {
            mismatch_samples.push_back(case_label + " finite mismatch: hmmer=" + std::to_string(hmmer_score) +
                                       " ours=" + std::to_string(our_score) +
                                       " abs_error=" + std::to_string(abs_error));
          }
        }
      } else {
        const bool same_inf = both_same_infinity(hmmer_score, our_score);
        if (same_inf) {
          ++nonfinite_equal_pairs;
          ++within_tolerance_pairs;
        } else {
          ++mismatched_pairs;
          ++nonfinite_mismatch_pairs;
          if (mismatch_samples.size() < 5) {
            mismatch_samples.push_back(case_label + " non-finite mismatch: hmmer=" + std::to_string(hmmer_score) +
                                       " ours=" + std::to_string(our_score));
          }
        }
      }

      ++processed_pairs;

      msv_profilectx_destroy(profile_ctx);
      msv_seqrecord_destroy(seq_record);
    }

    msv_seq_reader_close(seq_reader);
    msv_hmmrecord_destroy(hmm_record);
  }

  msv_hmm_reader_close(hmm_reader);

  if (finite_compare_pairs > 0) {
    mean_abs_error = error_sum / static_cast<float>(finite_compare_pairs);
  } else {
    mean_abs_error = std::numeric_limits<float>::quiet_NaN();
  }

  std::cout << "[integration] generic p7_GMSV parity summary"
            << " | hmms=" << processed_hmms
            << " pairs=" << processed_pairs
            << " within_tol=" << within_tolerance_pairs
            << " finite_pairs=" << finite_compare_pairs
            << " finite_mismatches=" << finite_mismatch_pairs
            << " nonfinite_equal=" << nonfinite_equal_pairs
            << " nonfinite_mismatches=" << nonfinite_mismatch_pairs
            << " mismatches=" << mismatched_pairs
            << " max_abs_error=" << max_abs_error
            << " mean_abs_error=" << mean_abs_error;
  if (!worst_case_label.empty()) {
    std::cout << " worst_case=" << worst_case_label;
  }
  std::cout << "\n";

  for (const std::string& sample : mismatch_samples) {
    std::cout << "[integration] sample mismatch: " << sample << "\n";
  }

  ASSERT_GT(processed_hmms, 0) << "No HMM records were processed";
  ASSERT_GT(processed_pairs, 0) << "No HMM/sequence pairs were processed";
  EXPECT_EQ(0, mismatched_pairs) << "Detected mismatches vs generic p7_GMSV; see summary above";
}

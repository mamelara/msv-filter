#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "aa_alphabet.hpp"
#include "hmmer_c_bridge.h"
#include "msv_bridge_adapter.hpp"
#include "msv_filter.hpp"

namespace {

struct RunConfig {
  std::string hmm_path;
  std::string fasta_path;
  int max_hmms;
  int max_seqs;
  float expected_hit_count;
};

struct RunStats {
  int processed_hmms = 0;
  int processed_pairs = 0;
  int bridge_failures = 0;
  int filter_failures = 0;
  int finite_scores = 0;
  int nonfinite_scores = 0;
  float sum_scores = 0.0f;
  float min_score = std::numeric_limits<float>::infinity();
  float max_score = -std::numeric_limits<float>::infinity();
  std::vector<std::string> failure_samples;
};

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

std::string label_for_case(const char* hmm_name, const char* seq_name, int hmm_idx, int seq_idx) {
  return "hmm#" + std::to_string(hmm_idx) + "(" + hmm_name + ")" +
         " vs seq#" + std::to_string(seq_idx) + "(" + seq_name + ")";
}

void add_failure_sample(RunStats* stats, const std::string& message) {
  if (stats->failure_samples.size() < 5) {
    stats->failure_samples.push_back(message);
  }
}

RunStats run_msv_filter_scan(const RunConfig& config) {
  RunStats stats;

  MsvHmmReader* hmm_reader = nullptr;
  const int open_status = msv_hmm_reader_open(config.hmm_path.c_str(), &hmm_reader);
  if (open_status != MSV_BRIDGE_OK || hmm_reader == nullptr) {
    ++stats.bridge_failures;
    add_failure_sample(&stats, "failed to open HMM path: " + config.hmm_path);
    return stats;
  }

  for (int hmm_idx = 0; hmm_idx < config.max_hmms; ++hmm_idx) {
    MsvHmmRecord* hmm_record = nullptr;
    const int next_hmm_status = msv_hmm_reader_next(hmm_reader, &hmm_record);
    if (next_hmm_status == MSV_BRIDGE_EOF) {
      break;
    }
    if (next_hmm_status != MSV_BRIDGE_OK || hmm_record == nullptr) {
      ++stats.bridge_failures;
      add_failure_sample(&stats, "failed reading HMM record: " + std::string(msv_hmm_reader_last_error(hmm_reader)));
      continue;
    }

    ++stats.processed_hmms;

    MsvSeqReader* seq_reader = nullptr;
    const int seq_open_status = msv_seq_reader_open(config.fasta_path.c_str(), hmm_record, &seq_reader);
    if (seq_open_status != MSV_BRIDGE_OK || seq_reader == nullptr) {
      ++stats.bridge_failures;
      add_failure_sample(&stats, "failed opening FASTA path: " + config.fasta_path);
      msv_hmmrecord_destroy(hmm_record);
      continue;
    }

    for (int seq_idx = 0; seq_idx < config.max_seqs; ++seq_idx) {
      MsvSeqRecord* seq_record = nullptr;
      const int next_seq_status = msv_seq_reader_next(seq_reader, &seq_record);
      if (next_seq_status == MSV_BRIDGE_EOF) {
        break;
      }
      if (next_seq_status != MSV_BRIDGE_OK || seq_record == nullptr) {
        ++stats.bridge_failures;
        add_failure_sample(&stats, "failed reading sequence: " + std::string(msv_seq_reader_last_error(seq_reader)));
        continue;
      }

      const std::string case_label = label_for_case(msv_hmmrecord_name(hmm_record),
                                                    msv_seqrecord_name(seq_record),
                                                    hmm_idx,
                                                    seq_idx);

      MsvProfileCtx* profile_ctx = nullptr;
      const int seq_length = static_cast<int>(msv_seqrecord_length(seq_record));
      const int ctx_status = msv_profilectx_create(hmm_record, seq_length, MSV_PROFILE_MODE_LOCAL, &profile_ctx);
      if (ctx_status != MSV_BRIDGE_OK || profile_ctx == nullptr) {
        ++stats.bridge_failures;
        add_failure_sample(&stats, "failed creating profile context for " + case_label);
        msv_seqrecord_destroy(seq_record);
        continue;
      }

      try {
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
            config.expected_hit_count,
            &our_score);

        if (our_status != 0) {
          ++stats.filter_failures;
          add_failure_sample(&stats, case_label + " status mismatch: ours=" + std::to_string(our_status));
        }

        if (std::isfinite(our_score)) {
          stats.sum_scores += our_score;
          ++stats.finite_scores;
          stats.min_score = std::min(stats.min_score, our_score);
          stats.max_score = std::max(stats.max_score, our_score);
        } else {
          ++stats.nonfinite_scores;
        }
      } catch (const std::exception& ex) {
        ++stats.bridge_failures;
        add_failure_sample(&stats, case_label + " adapter failure: " + ex.what());
      }

      ++stats.processed_pairs;

      msv_profilectx_destroy(profile_ctx);
      msv_seqrecord_destroy(seq_record);
    }

    msv_seq_reader_close(seq_reader);
    msv_hmmrecord_destroy(hmm_record);
  }

  msv_hmm_reader_close(hmm_reader);
  return stats;
}

}  // namespace

int main() {
  RunConfig config;
  config.hmm_path = get_env_or_default("MSV_ITEST_HMM_PATH", "databases/panther.100.hmm");
  config.fasta_path = get_env_or_default("MSV_ITEST_FASTA_PATH", "inputs/Arabidopsis_thaliana.100.pep.fa");
  config.max_hmms = get_env_int_or_default("MSV_ITEST_MAX_HMMS", 3);
  config.max_seqs = get_env_int_or_default("MSV_ITEST_MAX_SEQS", 25);
  config.expected_hit_count = get_env_float_or_default("MSV_ITEST_NU", 2.0f);

  if (config.max_hmms <= 0 || config.max_seqs <= 0) {
    std::cerr << "MSV_ITEST_MAX_HMMS and MSV_ITEST_MAX_SEQS must be > 0\n";
    return EXIT_FAILURE;
  }

  const RunStats stats = run_msv_filter_scan(config);

  std::cout << "[scan] msv_filter summary"
            << " | hmms=" << stats.processed_hmms
            << " pairs=" << stats.processed_pairs
            << " bridge_failures=" << stats.bridge_failures
            << " filter_failures=" << stats.filter_failures
            << " finite_scores=" << stats.finite_scores
            << " nonfinite_scores=" << stats.nonfinite_scores
            << " sum_scores=" << stats.sum_scores;

  if (stats.finite_scores > 0) {
    std::cout << " min_score=" << stats.min_score
              << " max_score=" << stats.max_score;
  }

  std::cout << "\n";

  for (const std::string& sample : stats.failure_samples) {
    std::cout << "[scan] sample failure: " << sample << "\n";
  }

  return EXIT_SUCCESS;
}

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

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
  int repeats;
  bool print_matrix;
  float expected_hit_count;
};

using HmmRecordPtr = std::unique_ptr<MsvHmmRecord, decltype(&msv_hmmrecord_destroy)>;
using ProfileCtxPtr = std::unique_ptr<MsvProfileCtx, decltype(&msv_profilectx_destroy)>;
using HmmReaderPtr = std::unique_ptr<MsvHmmReader, decltype(&msv_hmm_reader_close)>;

struct LoadedSequence {
  std::string name;
  int length = 0;
  std::vector<DigitalResidue> residues;
};

struct LoadedHmms {
  HmmReaderPtr reader{nullptr, &msv_hmm_reader_close};
  std::vector<HmmRecordPtr> records;
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

LoadedHmms load_hmms(const RunConfig& config) {
  MsvHmmReader* raw_reader = nullptr;
  const int open_status = msv_hmm_reader_open(config.hmm_path.c_str(), &raw_reader);
  if (open_status != MSV_BRIDGE_OK || raw_reader == nullptr) {
    throw std::runtime_error("failed to open HMM path: " + config.hmm_path);
  }

  LoadedHmms hmms;
  hmms.reader.reset(raw_reader);
  hmms.records.reserve(static_cast<size_t>(config.max_hmms));

  for (int hmm_idx = 0; hmm_idx < config.max_hmms; ++hmm_idx) {
    MsvHmmRecord* raw_hmm = nullptr;
    const int next_status = msv_hmm_reader_next(hmms.reader.get(), &raw_hmm);
    if (next_status == MSV_BRIDGE_EOF) {
      break;
    }
    if (next_status != MSV_BRIDGE_OK || raw_hmm == nullptr) {
      throw std::runtime_error("failed reading HMM record: " + std::string(msv_hmm_reader_last_error(hmms.reader.get())));
    }

    hmms.records.emplace_back(raw_hmm, &msv_hmmrecord_destroy);
  }

  if (hmms.records.empty()) {
    throw std::runtime_error("no HMM records were loaded from: " + config.hmm_path);
  }

  return hmms;
}

std::vector<LoadedSequence> load_sequences(const RunConfig& config, const MsvHmmRecord* hmm_record) {
  using SeqReaderPtr = std::unique_ptr<MsvSeqReader, decltype(&msv_seq_reader_close)>;
  using SeqRecordPtr = std::unique_ptr<MsvSeqRecord, decltype(&msv_seqrecord_destroy)>;

  MsvSeqReader* raw_reader = nullptr;
  const int open_status = msv_seq_reader_open(config.fasta_path.c_str(), hmm_record, &raw_reader);
  if (open_status != MSV_BRIDGE_OK || raw_reader == nullptr) {
    throw std::runtime_error("failed to open FASTA path: " + config.fasta_path);
  }

  SeqReaderPtr reader(raw_reader, &msv_seq_reader_close);
  std::vector<LoadedSequence> sequences;
  sequences.reserve(static_cast<size_t>(config.max_seqs));

  for (int seq_idx = 0; seq_idx < config.max_seqs; ++seq_idx) {
    MsvSeqRecord* raw_seq = nullptr;
    const int next_status = msv_seq_reader_next(reader.get(), &raw_seq);
    if (next_status == MSV_BRIDGE_EOF) {
      break;
    }
    if (next_status != MSV_BRIDGE_OK || raw_seq == nullptr) {
      throw std::runtime_error("failed reading sequence: " + std::string(msv_seq_reader_last_error(reader.get())));
    }

    SeqRecordPtr seq_record(raw_seq, &msv_seqrecord_destroy);
    sequences.push_back(LoadedSequence{
        msv_seqrecord_name(seq_record.get()),
        static_cast<int>(msv_seqrecord_length(seq_record.get())),
        make_digital_sequence_from_bridge(seq_record.get())});
  }

  if (sequences.empty()) {
    throw std::runtime_error("no sequence records were loaded from: " + config.fasta_path);
  }

  return sequences;
}

std::vector<std::vector<float>> run_msv_filter_scan(const RunConfig& config) {
  const LoadedHmms hmms = load_hmms(config);
  const std::vector<LoadedSequence> sequences = load_sequences(config, hmms.records.front().get());

  AminoAcidAlphabet alphabet;
  std::vector<std::vector<float>> score_matrix(
      hmms.records.size(), std::vector<float>(sequences.size(), 0.0f));


  for (size_t hmm_idx = 0; hmm_idx < hmms.records.size(); ++hmm_idx) {
    const MsvHmmRecord* hmm_record = hmms.records[hmm_idx].get();

    for (size_t seq_idx = 0; seq_idx < sequences.size(); ++seq_idx) {
      const LoadedSequence& sequence = sequences[seq_idx];

      MsvProfileCtx* raw_profile_ctx = nullptr;
      const int ctx_status = msv_profilectx_create(
          hmm_record, sequence.length, MSV_PROFILE_MODE_LOCAL, &raw_profile_ctx);
      if (ctx_status != MSV_BRIDGE_OK || raw_profile_ctx == nullptr) {
        throw std::runtime_error(
            "failed creating profile context for " +
            std::string(msv_hmmrecord_name(hmm_record)) + " vs " + sequence.name);
      }

      ProfileCtxPtr profile_ctx(raw_profile_ctx, &msv_profilectx_destroy);
      HMMProfile cpp_profile = make_profile_from_bridge(profile_ctx.get(), alphabet);
      DPMatrix cpp_dp(msv_hmmrecord_model_length(hmm_record), sequence.length);

      float score = 0.0f;
      const int status = msv_filter(
          sequence.residues.data(),
          sequence.length,
          cpp_profile,
          cpp_dp,
          config.expected_hit_count,
          &score);
      if (status != 0) {
        throw std::runtime_error(
            "msv_filter failed for " + std::string(msv_hmmrecord_name(hmm_record)) +
            " vs " + sequence.name);
      }

      score_matrix[hmm_idx][seq_idx] = score;
    }
  }

  return score_matrix;
}

}  // namespace

int main() {
  RunConfig config;
  config.hmm_path = get_env_or_default("MSV_ITEST_HMM_PATH", "databases/panther.100.hmm");
  config.fasta_path = get_env_or_default("MSV_ITEST_FASTA_PATH", "inputs/Arabidopsis_thaliana.100.pep.fa");
  config.max_hmms = get_env_int_or_default("MSV_ITEST_MAX_HMMS", 3);
  config.max_seqs = get_env_int_or_default("MSV_ITEST_MAX_SEQS", 25);
  config.repeats = get_env_int_or_default("MSV_SCAN_REPEATS", 1);
  config.print_matrix = get_env_int_or_default("MSV_SCAN_PRINT_MATRIX", 0) != 0;
  config.expected_hit_count = get_env_float_or_default("MSV_ITEST_NU", 2.0f);

  if (config.max_hmms <= 0 || config.max_seqs <= 0 || config.repeats <= 0) {
    std::cerr << "MSV_ITEST_MAX_HMMS, MSV_ITEST_MAX_SEQS, and MSV_SCAN_REPEATS must be > 0\n";
    return EXIT_FAILURE;
  }

  try {
    std::vector<std::vector<float>> score_matrix;
    const double start_time = omp_get_wtime();
    for (int repeat_idx = 0; repeat_idx < config.repeats; ++repeat_idx) {
      score_matrix = run_msv_filter_scan(config);
    }
    const double elapsed_seconds = omp_get_wtime() - start_time;

    std::cout << "[scan] score matrix rows=" << score_matrix.size();
    if (!score_matrix.empty()) {
      std::cout << " cols=" << score_matrix.front().size();
    }
    std::cout << " repeats=" << config.repeats
              << " elapsed_seconds=" << elapsed_seconds
              << " average_seconds=" << (elapsed_seconds / static_cast<double>(config.repeats))
              << "\n";

    if (config.print_matrix) {
      for (const std::vector<float>& row : score_matrix) {
        for (size_t col = 0; col < row.size(); ++col) {
          if (col > 0) {
            std::cout << ' ';
          }
          std::cout << row[col];
        }
        std::cout << "\n";
      }
    }
  } catch (const std::exception& ex) {
    std::cerr << "msv_filter_scan failed: " << ex.what() << "\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

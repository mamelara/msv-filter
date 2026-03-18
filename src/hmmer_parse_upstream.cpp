#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "aa_alphabet.hpp"
#include "hmmer_c_bridge.h"
#include "msv_bridge_adapter.hpp"

namespace {

void print_usage(const char* program_name) {
  std::cerr << "Usage: " << program_name << " <hmm-file> <fasta-file>\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  const char* hmm_path = argv[1];
  const char* fasta_path = argv[2];

  MsvHmmReader* hmm_reader = nullptr;
  MsvHmmRecord* hmm_record = nullptr;
  MsvSeqReader* seq_reader = nullptr;
  MsvSeqRecord* seq_record = nullptr;
  MsvProfileCtx* profile_ctx = nullptr;

  if (msv_hmm_reader_open(hmm_path, &hmm_reader) != MSV_BRIDGE_OK) {
    std::cerr << "Failed to open HMM file: " << hmm_path << "\n";
    return EXIT_FAILURE;
  }

  if (msv_hmm_reader_next(hmm_reader, &hmm_record) != MSV_BRIDGE_OK) {
    std::cerr << "Failed to read first HMM: " << msv_hmm_reader_last_error(hmm_reader) << "\n";
    msv_hmm_reader_close(hmm_reader);
    return EXIT_FAILURE;
  }

  if (msv_seq_reader_open(fasta_path, hmm_record, &seq_reader) != MSV_BRIDGE_OK) {
    std::cerr << "Failed to open FASTA: " << fasta_path << "\n";
    msv_hmmrecord_destroy(hmm_record);
    msv_hmm_reader_close(hmm_reader);
    return EXIT_FAILURE;
  }

  if (msv_seq_reader_next(seq_reader, &seq_record) != MSV_BRIDGE_OK) {
    std::cerr << "Failed to read first sequence: " << msv_seq_reader_last_error(seq_reader) << "\n";
    msv_seq_reader_close(seq_reader);
    msv_hmmrecord_destroy(hmm_record);
    msv_hmm_reader_close(hmm_reader);
    return EXIT_FAILURE;
  }

  const int sequence_length = static_cast<int>(msv_seqrecord_length(seq_record));
  if (msv_profilectx_create(hmm_record, sequence_length, MSV_PROFILE_MODE_LOCAL, &profile_ctx) != MSV_BRIDGE_OK) {
    std::cerr << "Failed to configure profile context\n";
    msv_seqrecord_destroy(seq_record);
    msv_seq_reader_close(seq_reader);
    msv_hmmrecord_destroy(hmm_record);
    msv_hmm_reader_close(hmm_reader);
    return EXIT_FAILURE;
  }

  try {
    AminoAcidAlphabet cpp_alphabet;
    auto digital_sequence = make_digital_sequence_from_bridge(seq_record);
    HMMProfile cpp_profile = make_profile_from_bridge(profile_ctx, cpp_alphabet);
    DPMatrix cpp_dp = make_dp_matrix_from_bridge(seq_record, profile_ctx);

    float baseline_score = 0.0f;
    if (msv_profilectx_score_gmsv(profile_ctx, seq_record, 2.0f, &baseline_score) != MSV_BRIDGE_OK) {
      throw std::runtime_error(msv_profilectx_last_error(profile_ctx));
    }

    std::cout << "Parsed HMM: " << msv_hmmrecord_name(hmm_record)
              << " (M=" << cpp_profile.model_length << ")\n";
    std::cout << "Parsed sequence: " << msv_seqrecord_name(seq_record)
              << " (L=" << sequence_length << ")\n";
    std::cout << "C++ residues buffer size: " << digital_sequence.size() << "\n";
    std::cout << "C++ DP matrix rows: " << cpp_dp.sequence_length + 1 << "\n";
    std::cout << "Baseline p7_GMSV score: " << baseline_score << "\n";
  } catch (const std::exception& ex) {
    std::cerr << "Bridge adapter failure: " << ex.what() << "\n";
    msv_profilectx_destroy(profile_ctx);
    msv_seqrecord_destroy(seq_record);
    msv_seq_reader_close(seq_reader);
    msv_hmmrecord_destroy(hmm_record);
    msv_hmm_reader_close(hmm_reader);
    return EXIT_FAILURE;
  }

  msv_profilectx_destroy(profile_ctx);
  msv_seqrecord_destroy(seq_record);
  msv_seq_reader_close(seq_reader);
  msv_hmmrecord_destroy(hmm_record);
  msv_hmm_reader_close(hmm_reader);
  return EXIT_SUCCESS;
}

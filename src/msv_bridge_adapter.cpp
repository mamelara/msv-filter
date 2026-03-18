#include "msv_bridge_adapter.hpp"

#include <algorithm>
#include <stdexcept>

std::vector<DigitalResidue> make_digital_sequence_from_bridge(const MsvSeqRecord* seq_record) {
  if (seq_record == nullptr) {
    throw std::invalid_argument("sequence record is null");
  }

  const int64_t sequence_length = msv_seqrecord_length(seq_record);
  if (sequence_length < 0) {
    throw std::invalid_argument("bridge returned an invalid sequence length");
  }

  std::vector<DigitalResidue> residues(static_cast<size_t>(sequence_length + 2));
  for (int i = 0; i <= sequence_length + 1; ++i) {
    residues[static_cast<size_t>(i)] = static_cast<DigitalResidue>(msv_seqrecord_residue_code(seq_record, i));
  }

  return residues;
}

HMMProfile make_profile_from_bridge(const MsvProfileCtx* profile_ctx, const AminoAcidAlphabet& alphabet) {
  if (profile_ctx == nullptr) {
    throw std::invalid_argument("profile context is null");
  }

  const int model_length = msv_profilectx_model_length(profile_ctx);
  const int target_length = msv_profilectx_target_length(profile_ctx);
  const int source_alphabet_size = msv_profilectx_alphabet_size(profile_ctx);

  if (model_length < 0 || target_length < 0 || source_alphabet_size < 0) {
    throw std::invalid_argument("bridge returned invalid profile metadata");
  }
  if (source_alphabet_size > alphabet.Kp) {
    throw std::invalid_argument("C++ alphabet does not cover all bridge residue codes");
  }

  HMMProfile profile(model_length, &alphabet);
  profile.model_length = model_length;
  profile.sequence_length = target_length;
  profile.mode = msv_profilectx_mode(profile_ctx);

  for (int k = 1; k <= model_length; ++k) {
    for (int residue = 0; residue < source_alphabet_size; ++residue) {
      profile.match_score(k, residue) = msv_profilectx_match_score(profile_ctx, k, residue);
    }
  }

  return profile;
}

DPMatrix make_dp_matrix_from_bridge(const MsvSeqRecord* seq_record, const MsvProfileCtx* profile_ctx) {
  if (seq_record == nullptr || profile_ctx == nullptr) {
    throw std::invalid_argument("sequence record and profile context are required");
  }

  const int64_t sequence_length = msv_seqrecord_length(seq_record);
  const int model_length = msv_profilectx_model_length(profile_ctx);
  if (sequence_length < 0 || model_length < 0) {
    throw std::invalid_argument("bridge returned invalid model/sequence lengths");
  }

  return DPMatrix(model_length, static_cast<int>(sequence_length));
}

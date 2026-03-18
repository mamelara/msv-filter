#ifndef MSV_FILTER_MSV_BRIDGE_ADAPTER_HPP
#define MSV_FILTER_MSV_BRIDGE_ADAPTER_HPP

#include <vector>

#include "aa_alphabet.hpp"
#include "dp_matrix.hpp"
#include "hmmer_c_bridge.h"
#include "hmmer_types.hpp"
#include "profile.hpp"

std::vector<DigitalResidue> make_digital_sequence_from_bridge(const MsvSeqRecord* seq_record);
HMMProfile make_profile_from_bridge(const MsvProfileCtx* profile_ctx, const AminoAcidAlphabet& alphabet);
DPMatrix make_dp_matrix_from_bridge(const MsvSeqRecord* seq_record, const MsvProfileCtx* profile_ctx);

#endif

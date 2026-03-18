#ifndef MSV_FILTER_HMMER_C_BRIDGE_H
#define MSV_FILTER_HMMER_C_BRIDGE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MsvHmmReader MsvHmmReader;
typedef struct MsvHmmRecord MsvHmmRecord;
typedef struct MsvSeqReader MsvSeqReader;
typedef struct MsvSeqRecord MsvSeqRecord;
typedef struct MsvProfileCtx MsvProfileCtx;

typedef enum {
  MSV_BRIDGE_OK = 0,
  MSV_BRIDGE_EOF = 1,
  MSV_BRIDGE_ERR = 2,
  MSV_BRIDGE_INVALID = 3
} MsvBridgeStatus;

typedef enum {
  MSV_PROFILE_MODE_LOCAL = 1,
  MSV_PROFILE_MODE_GLOCAL = 2,
  MSV_PROFILE_MODE_UNILOCAL = 3,
  MSV_PROFILE_MODE_UNIGLOCAL = 4
} MsvProfileMode;

int msv_hmm_reader_open(const char* hmm_path, MsvHmmReader** out_reader);
void msv_hmm_reader_close(MsvHmmReader* reader);
int msv_hmm_reader_next(MsvHmmReader* reader, MsvHmmRecord** out_record);
const char* msv_hmm_reader_last_error(const MsvHmmReader* reader);

void msv_hmmrecord_destroy(MsvHmmRecord* record);
const char* msv_hmmrecord_name(const MsvHmmRecord* record);
int msv_hmmrecord_model_length(const MsvHmmRecord* record);
int msv_hmmrecord_alphabet_size(const MsvHmmRecord* record);

int msv_seq_reader_open(const char* seq_path, const MsvHmmRecord* hmm_record, MsvSeqReader** out_reader);
void msv_seq_reader_close(MsvSeqReader* reader);
int msv_seq_reader_next(MsvSeqReader* reader, MsvSeqRecord** out_record);
const char* msv_seq_reader_last_error(const MsvSeqReader* reader);

void msv_seqrecord_destroy(MsvSeqRecord* record);
const char* msv_seqrecord_name(const MsvSeqRecord* record);
int64_t msv_seqrecord_length(const MsvSeqRecord* record);
uint8_t msv_seqrecord_residue_code(const MsvSeqRecord* record, int position);

int msv_profilectx_create(const MsvHmmRecord* hmm_record,
                          int target_length,
                          int profile_mode,
                          MsvProfileCtx** out_profile_ctx);
void msv_profilectx_destroy(MsvProfileCtx* profile_ctx);
const char* msv_profilectx_last_error(const MsvProfileCtx* profile_ctx);
int msv_profilectx_model_length(const MsvProfileCtx* profile_ctx);
int msv_profilectx_alphabet_size(const MsvProfileCtx* profile_ctx);
int msv_profilectx_target_length(const MsvProfileCtx* profile_ctx);
int msv_profilectx_mode(const MsvProfileCtx* profile_ctx);
float msv_profilectx_match_score(const MsvProfileCtx* profile_ctx, int model_position, int residue_code);
int msv_profilectx_score_gmsv(MsvProfileCtx* profile_ctx,
                              const MsvSeqRecord* seq_record,
                              float expected_hit_count,
                              float* out_score);

#ifdef __cplusplus
}
#endif

#endif

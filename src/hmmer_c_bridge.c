#include "hmmer_c_bridge.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_sq.h>
#include <esl_sqio.h>
#include <hmmer.h>

enum { kErrorBufferSize = 256 };

struct MsvHmmReader {
  P7_HMMFILE* hfp;
  ESL_ALPHABET* abc;
  char last_error[kErrorBufferSize];
};

struct MsvHmmRecord {
  P7_HMM* hmm;
  const ESL_ALPHABET* abc;
};

struct MsvSeqReader {
  ESL_SQFILE* sqfp;
  const ESL_ALPHABET* abc;
  char last_error[kErrorBufferSize];
};

struct MsvSeqRecord {
  ESL_SQ* sq;
};

struct MsvProfileCtx {
  P7_BG* bg;
  P7_PROFILE* gm;
  P7_GMX* gx;
  char last_error[kErrorBufferSize];
};

static void set_error(char* destination, const char* message) {
  if (destination == NULL) {
    return;
  }
  if (message == NULL) {
    destination[0] = '\0';
    return;
  }

  snprintf(destination, kErrorBufferSize, "%s", message);
}

int msv_hmm_reader_open(const char* hmm_path, MsvHmmReader** out_reader) {
  if (hmm_path == NULL || out_reader == NULL) {
    return MSV_BRIDGE_INVALID;
  }

  MsvHmmReader* reader = (MsvHmmReader*)calloc(1, sizeof(MsvHmmReader));
  if (reader == NULL) {
    return MSV_BRIDGE_ERR;
  }

  int status = p7_hmmfile_Open(hmm_path, NULL, &reader->hfp, reader->last_error);
  if (status != eslOK) {
    free(reader);
    return MSV_BRIDGE_ERR;
  }

  *out_reader = reader;
  return MSV_BRIDGE_OK;
}

void msv_hmm_reader_close(MsvHmmReader* reader) {
  if (reader == NULL) {
    return;
  }

  p7_hmmfile_Close(reader->hfp);
  esl_alphabet_Destroy(reader->abc);
  free(reader);
}

int msv_hmm_reader_next(MsvHmmReader* reader, MsvHmmRecord** out_record) {
  if (reader == NULL || out_record == NULL) {
    return MSV_BRIDGE_INVALID;
  }

  P7_HMM* hmm = NULL;
  int status = p7_hmmfile_Read(reader->hfp, &reader->abc, &hmm);
  if (status == eslEOF) {
    return MSV_BRIDGE_EOF;
  }
  if (status != eslOK) {
    set_error(reader->last_error, "Failed to read HMM record");
    return MSV_BRIDGE_ERR;
  }

  MsvHmmRecord* record = (MsvHmmRecord*)calloc(1, sizeof(MsvHmmRecord));
  if (record == NULL) {
    p7_hmm_Destroy(hmm);
    set_error(reader->last_error, "Failed to allocate HMM record wrapper");
    return MSV_BRIDGE_ERR;
  }

  record->hmm = hmm;
  record->abc = hmm->abc;

  *out_record = record;
  return MSV_BRIDGE_OK;
}

const char* msv_hmm_reader_last_error(const MsvHmmReader* reader) {
  if (reader == NULL) {
    return "null hmm reader";
  }
  return reader->last_error;
}

void msv_hmmrecord_destroy(MsvHmmRecord* record) {
  if (record == NULL) {
    return;
  }

  p7_hmm_Destroy(record->hmm);
  free(record);
}

const char* msv_hmmrecord_name(const MsvHmmRecord* record) {
  if (record == NULL || record->hmm == NULL || record->hmm->name == NULL) {
    return "";
  }
  return record->hmm->name;
}

int msv_hmmrecord_model_length(const MsvHmmRecord* record) {
  if (record == NULL || record->hmm == NULL) {
    return -1;
  }
  return record->hmm->M;
}

int msv_hmmrecord_alphabet_size(const MsvHmmRecord* record) {
  if (record == NULL || record->abc == NULL) {
    return -1;
  }
  return record->abc->Kp;
}

int msv_seq_reader_open(const char* seq_path, const MsvHmmRecord* hmm_record, MsvSeqReader** out_reader) {
  if (seq_path == NULL || hmm_record == NULL || hmm_record->abc == NULL || out_reader == NULL) {
    return MSV_BRIDGE_INVALID;
  }

  MsvSeqReader* reader = (MsvSeqReader*)calloc(1, sizeof(MsvSeqReader));
  if (reader == NULL) {
    return MSV_BRIDGE_ERR;
  }

  reader->abc = hmm_record->abc;
  int status = esl_sqfile_OpenDigital(reader->abc, seq_path, eslSQFILE_UNKNOWN, NULL, &reader->sqfp);
  if (status != eslOK) {
    set_error(reader->last_error, "Failed to open sequence file");
    free(reader);
    return MSV_BRIDGE_ERR;
  }

  *out_reader = reader;
  return MSV_BRIDGE_OK;
}

void msv_seq_reader_close(MsvSeqReader* reader) {
  if (reader == NULL) {
    return;
  }

  esl_sqfile_Close(reader->sqfp);
  free(reader);
}

int msv_seq_reader_next(MsvSeqReader* reader, MsvSeqRecord** out_record) {
  if (reader == NULL || out_record == NULL) {
    return MSV_BRIDGE_INVALID;
  }

  ESL_SQ* sq = esl_sq_CreateDigital(reader->abc);
  if (sq == NULL) {
    set_error(reader->last_error, "Failed to allocate sequence object");
    return MSV_BRIDGE_ERR;
  }

  int status = esl_sqio_Read(reader->sqfp, sq);
  if (status == eslEOF) {
    esl_sq_Destroy(sq);
    return MSV_BRIDGE_EOF;
  }
  if (status != eslOK) {
    esl_sq_Destroy(sq);
    set_error(reader->last_error, "Failed to parse sequence record");
    return MSV_BRIDGE_ERR;
  }

  MsvSeqRecord* record = (MsvSeqRecord*)calloc(1, sizeof(MsvSeqRecord));
  if (record == NULL) {
    esl_sq_Destroy(sq);
    set_error(reader->last_error, "Failed to allocate sequence wrapper");
    return MSV_BRIDGE_ERR;
  }

  record->sq = sq;
  *out_record = record;
  return MSV_BRIDGE_OK;
}

const char* msv_seq_reader_last_error(const MsvSeqReader* reader) {
  if (reader == NULL) {
    return "null sequence reader";
  }
  return reader->last_error;
}

void msv_seqrecord_destroy(MsvSeqRecord* record) {
  if (record == NULL) {
    return;
  }

  esl_sq_Destroy(record->sq);
  free(record);
}

const char* msv_seqrecord_name(const MsvSeqRecord* record) {
  if (record == NULL || record->sq == NULL || record->sq->name == NULL) {
    return "";
  }
  return record->sq->name;
}

int64_t msv_seqrecord_length(const MsvSeqRecord* record) {
  if (record == NULL || record->sq == NULL) {
    return -1;
  }
  return record->sq->n;
}

uint8_t msv_seqrecord_residue_code(const MsvSeqRecord* record, int position) {
  if (record == NULL || record->sq == NULL || record->sq->dsq == NULL) {
    return (uint8_t)eslDSQ_ILLEGAL;
  }
  if (position < 0 || position > (int)(record->sq->n + 1)) {
    return (uint8_t)eslDSQ_ILLEGAL;
  }
  return record->sq->dsq[position];
}

int msv_profilectx_create(const MsvHmmRecord* hmm_record,
                          int target_length,
                          int profile_mode,
                          MsvProfileCtx** out_profile_ctx) {
  if (hmm_record == NULL || hmm_record->hmm == NULL || out_profile_ctx == NULL || target_length < 0) {
    return MSV_BRIDGE_INVALID;
  }

  MsvProfileCtx* ctx = (MsvProfileCtx*)calloc(1, sizeof(MsvProfileCtx));
  if (ctx == NULL) {
    return MSV_BRIDGE_ERR;
  }

  ctx->bg = p7_bg_Create(hmm_record->hmm->abc);
  if (ctx->bg == NULL) {
    set_error(ctx->last_error, "Failed to allocate P7_BG");
    msv_profilectx_destroy(ctx);
    return MSV_BRIDGE_ERR;
  }

  ctx->gm = p7_profile_Create(hmm_record->hmm->M, hmm_record->hmm->abc);
  if (ctx->gm == NULL) {
    set_error(ctx->last_error, "Failed to allocate P7_PROFILE");
    msv_profilectx_destroy(ctx);
    return MSV_BRIDGE_ERR;
  }

  int status = p7_ProfileConfig(hmm_record->hmm, ctx->bg, ctx->gm, target_length, profile_mode);
  if (status != eslOK) {
    set_error(ctx->last_error, "p7_ProfileConfig failed");
    msv_profilectx_destroy(ctx);
    return MSV_BRIDGE_ERR;
  }

  ctx->gx = p7_gmx_Create(hmm_record->hmm->M, target_length);
  if (ctx->gx == NULL) {
    set_error(ctx->last_error, "Failed to allocate P7_GMX");
    msv_profilectx_destroy(ctx);
    return MSV_BRIDGE_ERR;
  }

  *out_profile_ctx = ctx;
  return MSV_BRIDGE_OK;
}

void msv_profilectx_destroy(MsvProfileCtx* profile_ctx) {
  if (profile_ctx == NULL) {
    return;
  }

  p7_gmx_Destroy(profile_ctx->gx);
  p7_profile_Destroy(profile_ctx->gm);
  p7_bg_Destroy(profile_ctx->bg);
  free(profile_ctx);
}

const char* msv_profilectx_last_error(const MsvProfileCtx* profile_ctx) {
  if (profile_ctx == NULL) {
    return "null profile context";
  }
  return profile_ctx->last_error;
}

int msv_profilectx_model_length(const MsvProfileCtx* profile_ctx) {
  if (profile_ctx == NULL || profile_ctx->gm == NULL) {
    return -1;
  }
  return profile_ctx->gm->M;
}

int msv_profilectx_alphabet_size(const MsvProfileCtx* profile_ctx) {
  if (profile_ctx == NULL || profile_ctx->gm == NULL || profile_ctx->gm->abc == NULL) {
    return -1;
  }
  return profile_ctx->gm->abc->Kp;
}

int msv_profilectx_target_length(const MsvProfileCtx* profile_ctx) {
  if (profile_ctx == NULL || profile_ctx->gm == NULL) {
    return -1;
  }
  return profile_ctx->gm->L;
}

int msv_profilectx_mode(const MsvProfileCtx* profile_ctx) {
  if (profile_ctx == NULL || profile_ctx->gm == NULL) {
    return -1;
  }
  return profile_ctx->gm->mode;
}

float msv_profilectx_match_score(const MsvProfileCtx* profile_ctx, int model_position, int residue_code) {
  if (profile_ctx == NULL || profile_ctx->gm == NULL || profile_ctx->gm->abc == NULL || model_position < 0 || residue_code < 0) {
    return -eslINFINITY;
  }
  if (model_position > profile_ctx->gm->M || residue_code >= profile_ctx->gm->abc->Kp) {
    return -eslINFINITY;
  }
  return P7P_MSC(profile_ctx->gm, model_position, residue_code);
}

int msv_profilectx_score_gmsv(MsvProfileCtx* profile_ctx,
                              const MsvSeqRecord* seq_record,
                              float expected_hit_count,
                              float* out_score) {
  if (profile_ctx == NULL || profile_ctx->gm == NULL || profile_ctx->gx == NULL || seq_record == NULL ||
      seq_record->sq == NULL || out_score == NULL) {
    return MSV_BRIDGE_INVALID;
  }

  int status = p7_GMSV(seq_record->sq->dsq,
                       (int)seq_record->sq->n,
                       profile_ctx->gm,
                       profile_ctx->gx,
                       expected_hit_count,
                       out_score);
  if (status != eslOK) {
    set_error(profile_ctx->last_error, "p7_GMSV failed");
    return MSV_BRIDGE_ERR;
  }

  return MSV_BRIDGE_OK;
}

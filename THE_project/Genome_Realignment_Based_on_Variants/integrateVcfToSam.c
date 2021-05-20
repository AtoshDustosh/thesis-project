#include "integrateVcfToSam.h"

static char aux_appended_type = 'Z';
static char aux_appended_tag[3] = {'X', 'V', '\0'};

// Extract #(extension_Xpart) more bases when extracting ref sequence
static const int extension_lpart = 0;
static const int extension_rpart = 0;

static int integration_sv_min_len = 0;  // minimal length for a SV
static int integration_sv_max_len = 0;  // maximal length for a SV
static int integration_strategy = 0;

void check_files(Options *opts) {
  if (getSamFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete for variants integration - lak input "
            "files *.sam/bam\n");
    exit(EXIT_FAILURE);
  }
  if (getVcfFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete for variants integration - lak input "
            "files *.vcf/bcf\n");
    exit(EXIT_FAILURE);
  }
  if (getFaFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete for variants integration - lak input "
            "files *.fa/fna/fasta\n");
    exit(EXIT_FAILURE);
  }
  if (getOutputFile(opts) == NULL) {
    fprintf(stderr,
            "Warning: output file not set. Set as default output file: "
            "%s\n",
            default_outputFile);
    setOutputFile(opts, default_outputFile);
  }
}

static inline int ifCanIntegrateAllele(RecVcf_bplus *rv, int alleleIdx,
                                       int startPos, int endPos) {
  // Ignore REF allele
  if (alleleIdx == 0) return 0;
  const char *var_alt = rv_allele(rv, alleleIdx);
  const char *var_ref = rv_allele(rv, 0);
  int length_alt = strlen(var_alt);
  int length_ref = strlen(var_ref);
  int64_t varStartPos = rv_pos(rv);
  int64_t varEndPos = varStartPos + length_ref - 1;
  switch (integration_strategy) {
    case _OPT_INTEGRATION_SNPONLY: {
      if (length_alt >= integration_sv_min_len ||
          length_ref >= integration_sv_min_len) {
        return 0;
      } else {
        // This is based on the assumption that varEndPos > varStartPoos
        return varStartPos <= endPos && varEndPos >= startPos;
      }
      break;
    }
    case _OPT_INTEGRATION_SVONLY: {
      if (length_alt < integration_sv_min_len &&
          length_ref < integration_sv_min_len) {
        return 0;
      } else {
        // This is based on the assumption that varEndPos > varStartPoos
        return varStartPos <= endPos && varEndPos >= startPos;
      }
      break;
    }
    case _OPT_INTEGRATION_ALL: {
      // This is based on the assumption that varEndPos > varStartPoos
      return varStartPos <= endPos && varEndPos >= startPos;
      break;
    }
    default: {
      fprintf(stderr, "Error: integration strategy error. \n");
      exit(EXIT_FAILURE);
    }
  }
}

static inline int count_integrated_allele(RecVcf_bplus *rv, int64_t startPos,
                                          int64_t endPos) {
  // This is based on the assumption that varEndPos > varStartPoos
  int64_t varStartPos = rv_pos(rv);
  int cnt_integrated_allele = 0;
  // Calculate the affected length
  for (int i = 0; i < rv_alleleCnt(rv); i++) {
    switch (bcf_get_variant_type(rv_object(rv), i)) {
      case VCF_REF: {
        break;
      }
      case VCF_SNP: {
        if (varStartPos >= startPos && varStartPos <= endPos)
          cnt_integrated_allele++;
        break;
      }
      case VCF_INDEL:
      case VCF_MNP:
      case VCF_OTHER: {
        const char *var_ref = rv_allele(rv, 0);
        const char *var_alt = rv_allele(rv, i);
        int length_ref = strlen(var_ref);
        int length_alt = strlen(var_alt);
        int64_t varEndPos = varStartPos + length_ref - 1;
        assert((var_alt[0] != '<') ||
               (fprintf(stderr,
                        "Error: variants with tags should be ignored when "
                        "loading.\n") < 0));
        switch (integration_strategy) {
          case _OPT_INTEGRATION_SNPONLY: {
            if (length_alt < integration_sv_min_len &&
                length_ref < integration_sv_min_len) {
              if (varEndPos >= startPos && varStartPos <= endPos)
                cnt_integrated_allele++;
            }
            break;
          }
          case _OPT_INTEGRATION_SVONLY: {
            if (length_alt >= integration_sv_min_len ||
                length_ref >= integration_sv_min_len) {
              if (varEndPos >= startPos && varEndPos <= endPos)
                cnt_integrated_allele++;
            }
            break;
          }
          case _OPT_INTEGRATION_ALL: {
            if (varEndPos >= startPos && varStartPos <= endPos)
              cnt_integrated_allele++;
            break;
          }
          default: {
            fprintf(stderr, "Error: integration strategy error. \n");
            exit(EXIT_FAILURE);
          }
        }
        break;
      }
      default: {  // Ignore other kinds of variants (VCF_BND, VCF_OVERLAP)
        fprintf(stderr,
                "Warning: variant ignored. varRef: %s, varAlt[%d]: %s\n",
                rv_allele(rv, 0), i, rv_allele(rv, i));
      }
    }
  }
  return cnt_integrated_allele;
}

static inline AlignResult *integration_integrate_lpart(
    Element_RecVcf *ervArray[], int ervCombi[], int alleleCombi[],
    int length_combi, int64_t lbound_var, int64_t lbound_M, RecSam *rec_rs,
    GenomeFa *gf, GenomeSam *gs, GenomeVcf_bplus *gv, samFile *file_output,
    int *ret_length_lpart_ref) {
  // printf("lpart integration selected rv: \n");
  // for (int i = 0; i < length_combi; i++) {
  //   RecVcf_bplus *rv = ervArray[ervCombi[i]]->rv;
  //   genomeVcf_bplus_printRec(gv, rv);
  // }

  // Set bounds of ref sequence and update bounds
  // printf("selected alleles: \n");
  int64_t rbound_ref = lbound_M - 1;  // 1-based, included
  int64_t lbound_ref = lbound_var;
  int length_seq_ref = 0;
  for (int i = 0; i < length_combi; i++) {
    RecVcf_bplus *rv = ervArray[ervCombi[i]]->rv;
    const char *allele_ref = rv_allele(rv, 0);
    const char *allele_alt = rv_allele(rv, alleleCombi[i]);
    // printf("allele_ref[%d]: %s, allele_alt[%d]: %s\n", i, allele_ref, i,
    //        allele_alt);
    int length_allele_ref = strlen(allele_ref);
    int length_allele_alt = strlen(allele_alt);
    if (length_allele_ref == length_allele_alt) {
      if (length_allele_ref == 1) {
        continue;
      } else {
        length_seq_ref += length_allele_ref;
      }
    } else if (length_allele_alt > length_allele_ref) {
      length_seq_ref += (length_allele_alt - length_allele_ref);
    } else {
      lbound_ref -= length_allele_ref - length_allele_alt;
    }
  }
  lbound_ref -= extension_lpart;
  lbound_ref = lbound_ref <= 0 ? 1 : lbound_ref;

  // printf("lbound ref (before): %" PRId64 ", (after): %" PRId64 "\n",
  // lbound_var,
  //        lbound_ref);

  length_seq_ref += rbound_ref - lbound_ref + 1;
  if (length_seq_ref == 0) {
    // printf(">>>>>>>>>>>>>>>>>>>>>>>>>> empty lpart ref occurred. \n");
    // printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
    return init_AlignResult();
  }

  // A buffer for constructing new ref sequence
  char buf_seq[length_seq_ref + 1];  // "+1" indicates '\0' at the end
  memset(buf_seq, 0, length_seq_ref + 1);

  // Get original ref sequence
  const char *rname_read = rsDataRname(gs, rec_rs);
  ChromFa *tmp_cf = getChromFromGenomeFabyName(rname_read, gf);
  char *original_seq_ref = getSeqFromChromFa(lbound_ref, rbound_ref, tmp_cf);
  // printf("ref seq got [%" PRId64 ",%" PRId64 "]: %s\n", lbound_ref,
  // rbound_ref,
  //        original_seq_ref);

  // Integrate selected variants with original ref sequence
  int idx_seq_ref_new = 0;
  int idx_seq_ref_old = 0;
  for (int i = 0; i < length_combi; i++) {
    RecVcf_bplus *rv = ervArray[ervCombi[i]]->rv;
    int64_t pos_allele_start = rv_pos(rv);
    const char *allele_ref = rv_allele(rv, 0);
    const char *allele_alt = rv_allele(rv, alleleCombi[i]);
    // printf("allele_ref[%d]: %s, allele_alt[%d]: %s\n", i, allele_ref, i,
    // allele_alt);
    int length_allele_ref = strlen(allele_ref);
    int length_allele_alt = strlen(allele_alt);
    int idx_allele_ref = 0;
    int idx_allele_alt = 0;
    // Synchronize the positions of alleles before integration, especially when
    // the allele is a DEL and the start pos is not within the ref region.
    // Ignore the bases in front of the ref region.
    if (pos_allele_start < lbound_ref + idx_seq_ref_old) {
      idx_allele_alt += lbound_ref + idx_seq_ref_old - pos_allele_start;
      idx_allele_ref += lbound_ref + idx_seq_ref_old - pos_allele_start;
    }
    if (pos_allele_start >= rbound_ref) {
      break;
    }
    // Copy bases between this allele and last integrated allele on old ref seq
    while (idx_seq_ref_old + lbound_ref < pos_allele_start) {
      buf_seq[idx_seq_ref_new++] = original_seq_ref[idx_seq_ref_old++];
    }
    // Copy remained bases in ALT field
    while (idx_allele_alt < length_allele_alt) {
      assert(idx_seq_ref_new < length_seq_ref);
      buf_seq[idx_seq_ref_new++] = allele_alt[idx_allele_alt++];
    }
    // Ignore remained bases in REF field
    while (idx_allele_ref < length_allele_ref) {
      idx_seq_ref_old++;
      idx_allele_ref++;
    }
  }
  // Pad the remained unintegrated bases on old ref seq
  while (idx_seq_ref_old + lbound_ref <= rbound_ref) {
    buf_seq[idx_seq_ref_new++] = original_seq_ref[idx_seq_ref_old++];
  }
  buf_seq[idx_seq_ref_new] = '\0';
  // printf("integrated ref seq (L): %s\n", buf_seq);

  // Extract lpart ref seq
  char *seq_read = rsDataSeq(rec_rs);
  uint32_t length_read = rsDataSeqLength(rec_rs);
  int pos_start_lpart = 1;  // 1-based, included (lbound of M area)
  int length_M_area = 0;
  int tmp_pos_start_lpart = 1;
  // ATTENTION!!! See comments in rpart process.
  for (int i = 0; i < rs_cigar_cnt(rec_rs); i++) {
    uint32_t tmp_length = rs_cigar_oplen(rec_rs, i);
    char cigar_opChar = rs_cigar_opChar(rec_rs, i);
    if (cigar_opChar == 'M') {
      if (tmp_length >= length_M_area) {
        length_M_area = tmp_length;
        pos_start_lpart = tmp_pos_start_lpart;
      }
    }
    if (cigar_opChar == 'D' || cigar_opChar == 'H') {
      // Do not move on read sequence in such case
    } else {
      // For cigar op 'MX=IS', move on read sequence
      tmp_pos_start_lpart += tmp_length;
    }
  }
  int begin_subSeq = 0;
  // lbound_M - 1 (move out of the M area) - 1 (array index)
  int end_subSeq = pos_start_lpart - 2;
  char *seq_read_lpart = subStr(seq_read, begin_subSeq, end_subSeq);
  // printf("lpart begin: %d, end: %d\n", begin_subSeq, end_subSeq);
  // printf("extracted lpart seq: %s\n", seq_read_lpart);
  // seq_read_lpart can be NULL, which indicates that there is bases in lpart.
  if (seq_read_lpart == NULL) {
    free(seq_read);
    // printf(">>>>>>>>>>>>>>>>>>>>>>>>>> empty lpart seq occurred. \n");
    return init_AlignResult();
  }

  // Reverse ref seq and read seq lpart before alignment
  char *seq_read_lpart_rev = revStr(seq_read_lpart, strlen(seq_read_lpart));
  *ret_length_lpart_ref = strlen(buf_seq);
  char *seq_ref_lpart_rev = revStr(buf_seq, *ret_length_lpart_ref);

  // Align tseq and qseq using ksw2 global alignment
  AlignResult *ar = init_AlignResult();
  align_ksw2(seq_ref_lpart_rev, strlen(seq_ref_lpart_rev), seq_read_lpart_rev,
             strlen(seq_read_lpart_rev), ar);

  free(seq_read_lpart_rev);
  free(seq_ref_lpart_rev);

  free(seq_read_lpart);
  free(seq_read);

  // printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
  return ar;
}

static inline AlignResult *integration_integrate_rpart(
    Element_RecVcf *ervArray[], int ervCombi[], int alleleCombi[],
    int length_combi, int64_t rbound_M, int64_t rbound_var, RecSam *rec_rs,
    GenomeFa *gf, GenomeSam *gs, GenomeVcf_bplus *gv, samFile *file_output) {
  // printf("rpart integration selected rv: \n");
  // for (int i = 0; i < length_combi; i++) {
  //   RecVcf_bplus *rv = ervArray[ervCombi[i]]->rv;
  //   genomeVcf_bplus_printRec(gv, rv);
  // }

  // Set bounds of ref sequence and update bounds
  // printf("selected alleles: \n");
  int64_t lbound_ref = rbound_M + 1;  // 1-based, included
  int64_t rbound_ref = rbound_var;    // 1-based, included
  int length_seq_ref = 0;
  for (int i = 0; i < length_combi; i++) {
    const char *allele_ref = rv_allele(ervArray[ervCombi[i]]->rv, 0);
    const char *allele_alt =
        rv_allele(ervArray[ervCombi[i]]->rv, alleleCombi[i]);
    // printf("allele_ref[%d]: %s, allele_alt[%d]: %s\n", i, allele_ref, i,
    //        allele_alt);
    int length_allele_ref = strlen(allele_ref);
    int length_allele_alt = strlen(allele_alt);
    if (length_allele_alt == length_allele_ref) {
      if (length_allele_ref == 1) {
        continue;
      } else {
        length_seq_ref += length_allele_ref;
      }
    } else if (length_allele_alt > length_allele_ref) {
      length_seq_ref += (length_allele_alt - length_allele_ref);
    } else {
      rbound_ref += (length_allele_ref - length_allele_alt);
    }
  }
  rbound_ref += extension_rpart;

  // printf("rbound ref (before): %" PRId64 ", (after): %" PRId64 "\n",
  // rbound_var,
  //        rbound_ref);

  length_seq_ref += rbound_ref - lbound_ref + 1;
  if (length_seq_ref == 0) {
    // printf(">>>>>>>>>>>>>>>>>>>>>>>>>> empty rpart ref occurred. \n");
    AlignResult *ar = init_AlignResult();
    int length_M_area = 0;
    int idx_cigar_M = 0;
    // The following codes are used to locate the index of M area in the
    // original cigar and check if there exists 'I' in the end of the cigar. Pad
    // 'I' to the end of the align result if it exists and length_seq_ref == 0.
    // If this step is ignored, in some cases, bug will occur.
    for (int i = 0; i < rs_cigar_cnt(rec_rs); i++) {
      uint32_t tmp_length = rs_cigar_oplen(rec_rs, i);
      char cigar_opChar = rs_cigar_opChar(rec_rs, i);
      if (cigar_opChar == 'M') {
        // ATTENTION! This ">=" must be synchronized with the one used to find
        // bounds_M on reference. So does in the lpart process.
        if (tmp_length >= length_M_area) {
          length_M_area = tmp_length;
          idx_cigar_M = i;
        }
      }
    }
    // This is designed for sam records with cigar like "????3I" only. These
    // records will have empty rpart ref in the process, and the "#I" at the end
    // of original will be ignored, which will result in inconsistence betweeen
    // merged cigar and read sequence.
    if (idx_cigar_M == rs_cigar_cnt(rec_rs) - 2 && idx_cigar_M >= 0) {
      uint32_t cigar_len = rs_cigar_oplen(rec_rs, idx_cigar_M + 1);
      char cigar_opChar = rs_cigar_opChar(rec_rs, idx_cigar_M + 1);
      if (cigar_opChar == 'I') {
        ar->cigarLen = (uint32_t *)malloc(sizeof(uint32_t));
        ar->cigarOp = (char *)malloc(sizeof(char));
        ar->cnt_cigar = 1;
        ar->cigarLen[0] = cigar_len;
        ar->cigarOp[0] = cigar_opChar;
      }
    }
    // printf("###############################################################\n");
    return ar;
  }

  // A buffer for constructing new ref sequence
  char buf_seq[length_seq_ref + 1];  // "+1" indicates '\0' at the end
  memset(buf_seq, 0, length_seq_ref + 1);

  // Get original ref sequence
  const char *rname_read = rsDataRname(gs, rec_rs);
  ChromFa *tmp_cf = getChromFromGenomeFabyName(rname_read, gf);
  char *original_seq_ref = getSeqFromChromFa(lbound_ref, rbound_ref, tmp_cf);
  // printf("ref seq got [%" PRId64 ",%" PRId64 "]: %s\n", lbound_ref,
  // rbound_ref,
  //        original_seq_ref);

  // Integrate selected variants with original ref sequence
  int idx_seq_ref_new = 0;
  int idx_seq_ref_old = 0;
  for (int i = 0; i < length_combi; i++) {
    RecVcf_bplus *rv = ervArray[ervCombi[i]]->rv;
    int64_t pos_allele_start = rv_pos(rv);
    const char *allele_ref = rv_allele(rv, 0);
    const char *allele_alt = rv_allele(rv, alleleCombi[i]);
    // printf("allele_ref[%d]: %s, allele_alt[%d]: %s\n", i, allele_ref, i,
    //        allele_alt);
    int length_allele_ref = strlen(allele_ref);
    int length_allele_alt = strlen(allele_alt);
    int idx_allele_ref = 0;
    int idx_allele_alt = 0;
    // Synchronize the positions of alleles before integration, especially when
    // the allele is a DEL and the start pos is not within the ref region.
    // Ignore the bases in front of ref region.
    if (pos_allele_start < lbound_ref + idx_seq_ref_old) {
      idx_allele_alt += lbound_ref + idx_seq_ref_old - pos_allele_start;
      idx_allele_ref += lbound_ref + idx_seq_ref_old - pos_allele_start;
    }
    // Copy bases between this allele and last integrated allele on old ref seq
    while (idx_seq_ref_old + lbound_ref < pos_allele_start) {
      buf_seq[idx_seq_ref_new++] = original_seq_ref[idx_seq_ref_old++];
    }
    // Copy remained bases in ALT field
    while (idx_allele_alt < length_allele_alt) {
      assert(idx_seq_ref_new < length_seq_ref);
      buf_seq[idx_seq_ref_new++] = allele_alt[idx_allele_alt++];
    }
    // Ignore remained bases in REF field
    while (idx_allele_ref < length_allele_ref) {
      idx_seq_ref_old++;
      idx_allele_ref++;
    }
  }
  // Pad the remained unintegrated bases on old ref seq
  while (idx_seq_ref_old + lbound_ref <= rbound_ref) {
    buf_seq[idx_seq_ref_new++] = original_seq_ref[idx_seq_ref_old++];
  }
  buf_seq[idx_seq_ref_new] = '\0';
  // printf("integrated ref seq (R): %s\n", buf_seq);

  // Extract rpart read seq
  int pos_start_rpart = 1;  // 1-based, included (lbound of M area)
  int length_M_area = 0;
  int tmp_pos_start_rpart = 1;
  // The following codes are similar to the area in method "integration_threads"
  // where is commented with "find longest 'M' area in cigar". But they are
  // different. That one is used to find bounds_M on the reference. This one is
  // used to find bounds_M on the read seq. Despite they have different
  // purposes, their logics should be the same.
  for (int i = 0; i < rs_cigar_cnt(rec_rs); i++) {
    uint32_t tmp_length = rs_cigar_oplen(rec_rs, i);
    char cigar_opChar = rs_cigar_opChar(rec_rs, i);
    if (cigar_opChar == 'M') {
      // ATTENTION! This ">=" must be synchronized with the one used to find
      // bounds_M on reference. So does in the lpart process.
      if (tmp_length >= length_M_area) {
        length_M_area = tmp_length;
        pos_start_rpart = tmp_pos_start_rpart;
      }
    }
    if (cigar_opChar == 'D' || cigar_opChar == 'H') {
      // Do not move on read sequence in such case
    } else {
      // For cigar op 'MX=IS', move on read sequence
      tmp_pos_start_rpart += tmp_length;
    }
  }
  // lbound_M + (length_M - 1) + 1 (move out of the M area) - 1 (array index)
  char *seq_read = rsDataSeq(rec_rs);
  uint32_t length_read = rsDataSeqLength(rec_rs);
  int begin_subSeq = pos_start_rpart + length_M_area - 1;
  int end_subSeq = length_read - 1;
  char *seq_read_rpart = subStr(seq_read, begin_subSeq, end_subSeq);
  // printf("rpart begin: %d, end: %d\n", begin_subSeq, end_subSeq);
  // printf("extracted rpart seq: %s\n", seq_read_rpart);
  if (seq_read_rpart == NULL) {
    free(seq_read);
    // printf(">>>>>>>>>>>>>>>>>>>>>>>>>> empty rpart seq occurred. \n");
    return init_AlignResult();
  }

  // Align tseq and qseq using ksw2 global alignment
  AlignResult *ar = init_AlignResult();
  align_ksw2(buf_seq, strlen(buf_seq), seq_read_rpart, strlen(seq_read_rpart),
             ar);

  free(seq_read_rpart);
  free(seq_read);

  // printf("###############################################################\n");
  return ar;
}

static inline void integration_integrate(
    Element_RecVcf *ervArray_lpart[], int ervCombi_lpart[],
    int alleleCombi_lpart[], int length_combi_lpart, int length_ervArray_lpart,
    Element_RecVcf *ervArray_rpart[], int ervCombi_rpart[],
    int alleleCombi_rpart[], int length_combi_rpart, int length_ervArray_rpart,
    int64_t lbound_var, int64_t rbound_var, int64_t lbound_M, int64_t rbound_M,
    RecSam *rec_rs, int64_t id_rec, GenomeFa *gf, GenomeSam *gs,
    GenomeVcf_bplus *gv, samFile *file_output) {
  // Integrate the left part
  int ret_length_lpart_ref = 0;
  AlignResult *ar_lpart = integration_integrate_lpart(
      ervArray_lpart, ervCombi_lpart, alleleCombi_lpart, length_combi_lpart,
      lbound_var, lbound_M, rec_rs, gf, gs, gv, file_output,
      &ret_length_lpart_ref);
  // Integrate the right part
  AlignResult *ar_rpart = integration_integrate_rpart(
      ervArray_rpart, ervCombi_rpart, alleleCombi_rpart, length_combi_rpart,
      rbound_M, rbound_var, rec_rs, gf, gs, gv, file_output);

  // Fix cigars: remove leftmost 'D' and rightmost 'D'
  // And calculate new POS for the alignment result
  int64_t new_pos = lbound_M - ret_length_lpart_ref;
  int cnt_cigar_lpart = ar_cigar_cnt(ar_lpart);
  int cnt_cigar_rpart = ar_cigar_cnt(ar_rpart);
  if (cnt_cigar_lpart > 0 && ar_cigarOp(ar_lpart, cnt_cigar_lpart - 1) == 'D') {
    new_pos += ar_cigarlen(ar_lpart, cnt_cigar_lpart - 1);
    ar_lpart->cnt_cigar--;
    cnt_cigar_lpart--;
  }
  if (cnt_cigar_rpart > 0 && ar_cigarOp(ar_rpart, cnt_cigar_rpart - 1) == 'D') {
    ar_rpart->cnt_cigar--;
    cnt_cigar_rpart--;
  }

  // Merge cigar
  uint32_t length_M = rbound_M - lbound_M + 1;
  // printf("length M area: %" PRIu32 "\n", length_M);
  int length_buf_cigar = cnt_cigar_lpart + cnt_cigar_rpart + 1;
  int idx_buf_cigar = 0;
  uint32_t buf_cigarLen[length_buf_cigar];
  char buf_cigarOp[length_buf_cigar];
  length_buf_cigar = 0;  // recount length of cigar
  for (int i = cnt_cigar_lpart - 1; i >= 0; i--) {
    buf_cigarLen[idx_buf_cigar] = ar_cigarlen(ar_lpart, i);
    buf_cigarOp[idx_buf_cigar] = ar_cigarOp(ar_lpart, i);
    idx_buf_cigar++;
    length_buf_cigar++;
  }
  if (idx_buf_cigar > 0) {
    if (buf_cigarOp[idx_buf_cigar - 1] == 'M') {
      buf_cigarLen[idx_buf_cigar - 1] += length_M;
    } else {
      buf_cigarOp[idx_buf_cigar] = 'M';
      buf_cigarLen[idx_buf_cigar] = length_M;
      idx_buf_cigar++;
      length_buf_cigar++;
    }
  } else {
    buf_cigarOp[idx_buf_cigar] = 'M';
    buf_cigarLen[idx_buf_cigar] = length_M;
    idx_buf_cigar++;
    length_buf_cigar++;
  }
  for (int i = 0; i < cnt_cigar_rpart; i++) {
    if (i == 0) {
      if (ar_cigarOp(ar_rpart, 0) == 'M') {
        buf_cigarLen[idx_buf_cigar - 1] += ar_cigarlen(ar_rpart, 0);
      } else {
        buf_cigarLen[idx_buf_cigar] = ar_cigarlen(ar_rpart, 0);
        buf_cigarOp[idx_buf_cigar] = ar_cigarOp(ar_rpart, 0);
        idx_buf_cigar++;
        length_buf_cigar++;
      }
    } else {
      buf_cigarLen[idx_buf_cigar] = ar_cigarlen(ar_rpart, i);
      buf_cigarOp[idx_buf_cigar] = ar_cigarOp(ar_rpart, i);
      idx_buf_cigar++;
      length_buf_cigar++;
    }
  }
  // "* 8": you can't put a cigar op with cigar len bigger than 10^8
  char buf_merged_cigar[length_buf_cigar * 8];
  memset(buf_merged_cigar, 0, length_buf_cigar * 8);
  for (int i = 0; i < length_buf_cigar; i++) {
    char tmp_cigar[10];
    sprintf(tmp_cigar, "%" PRIu32 "%c", buf_cigarLen[i], buf_cigarOp[i]);
    strcat(buf_merged_cigar, tmp_cigar);
  }

  // Print information collected for this integration task
  // printf("lbound var: %" PRId64 ", rbound var: %" PRId64 "\n", lbound_var,
  //        rbound_var);
  // printf("lbound_M: %" PRId64 ", rbound_M: %" PRId64 "\n", lbound_M,
  // rbound_M); printf("length combi lpart: %d, combi rpart: %d\n",
  // length_combi_lpart,
  //        length_combi_rpart);
  // printf("ervArray (L): \n");
  // for (int i = 0; i < length_ervArray_lpart; i++) {
  //   genomeVcf_bplus_printRec(gv, ervArray_lpart[i]->rv);
  // }
  // printf("\n");
  // printf("(rv, allele) selected (L): \n");
  // for (int i = 0; i < length_combi_lpart; i++) {
  //   printf("(%d,%d) ", ervCombi_lpart[i], alleleCombi_lpart[i]);
  // }
  // printf("\n");
  // printf("ervArray (R): \n");
  // for (int i = 0; i < length_ervArray_rpart; i++) {
  //   genomeVcf_bplus_printRec(gv, ervArray_rpart[i]->rv);
  // }
  // printf("\n");
  // printf("(rv, allele) selected (R): \n");
  // for (int i = 0; i < length_combi_rpart; i++) {
  //   printf("(%d,%d) ", ervCombi_rpart[i], alleleCombi_rpart[i]);
  // }
  // printf("\n");
  // printf("*****************************************************\n");

  // printf("lpart align");
  // print_AlignResult(ar_lpart);
  // printf("rpart align");
  // print_AlignResult(ar_rpart);

  // Write result into file
  // printSamRecord_brief(gs, rsData(rec_rs));
  // printf("merged cigar: %s\n", buf_merged_cigar);
  // printf("fixed POS: %" PRId64 "\n", new_pos);
  bam1_t *new_rec = bamSetPosCigarMapq(rsData(rec_rs), new_pos, 0,
                                       rsDataSeqLength(rec_rs) - 1,
                                       buf_merged_cigar, rsDataMapQ(rec_rs));
  // Add integrated variants' IDs of the lpart into the auxiliary fields
  int length_aux_data_lpart = 0;
  for (int i = 0; i < length_combi_lpart; i++) {
    RecVcf_bplus *rv = ervArray_lpart[ervCombi_lpart[i]]->rv;
    int64_t var_pos = rv_pos(rv);
    const char *str_id = rv_ID(rv);
    int idx_allele = alleleCombi_lpart[i];
    int length_str_id = strlen(str_id);
    char buf[length_str_id + 64];  // If there is more than 10^6 variants in
                                   // this combination, you will never see the
                                   // program finish.
    sprintf(buf, "%s;%" PRId64 ";%d ", str_id, var_pos, idx_allele);
    length_aux_data_lpart += strlen(buf);
  }
  char aux_data_lpart[length_aux_data_lpart + 1];
  memset(aux_data_lpart, 0, length_aux_data_lpart + 1);
  for (int i = 0; i < length_combi_lpart; i++) {
    RecVcf_bplus *rv = ervArray_lpart[ervCombi_lpart[i]]->rv;
    int64_t var_pos = rv_pos(rv);
    const char *str_id = rv_ID(rv);
    int idx_allele = alleleCombi_lpart[i];
    int length_str_id = strlen(str_id);
    char buf[length_str_id + 64];  // If there is more than 10^6 variants in
                                   // this combination, you will never see the
                                   // program finish.
    sprintf(buf, "%s;%" PRId64 ";%d ", str_id, var_pos, idx_allele);
    strcat(aux_data_lpart, buf);
  }
  int length_aux_data_rpart = 0;
  for (int i = 0; i < length_combi_rpart; i++) {
    RecVcf_bplus *rv = ervArray_rpart[ervCombi_rpart[i]]->rv;
    int64_t var_pos = rv_pos(rv);
    const char *str_id = rv_ID(rv);
    int idx_allele = alleleCombi_rpart[i];
    int length_str_id = strlen(str_id);
    char buf[length_str_id + 64];  // If there is more than 10^6 variants in
                                   // this combination, you will never see the
                                   // program finish.
    sprintf(buf, "%s;%" PRId64 ";%d ", str_id, var_pos, idx_allele);
    length_aux_data_rpart += strlen(buf);
  }
  char aux_data_rpart[length_aux_data_rpart + 1];
  memset(aux_data_rpart, 0, length_aux_data_rpart + 1);
  for (int i = 0; i < length_combi_rpart; i++) {
    RecVcf_bplus *rv = ervArray_rpart[ervCombi_rpart[i]]->rv;
    int64_t var_pos = rv_pos(rv);
    const char *str_id = rv_ID(rv);
    int idx_allele = alleleCombi_rpart[i];
    int length_str_id = strlen(str_id);
    char buf[length_str_id + 64];  // If there is more than 10^6 variants in
                                   // this combination, you will never see the
                                   // program finish.
    sprintf(buf, "%s;%" PRId64 ";%d ", str_id, var_pos, idx_allele);
    strcat(aux_data_rpart, buf);
  }
  int length_aux_data = length_aux_data_lpart + length_aux_data_rpart;
  char aux_data[length_aux_data + 1];
  memset(aux_data, 0, length_aux_data + 1);
  strcat(aux_data, aux_data_lpart);
  strcat(aux_data, aux_data_rpart);
  bam_aux_append(new_rec, aux_appended_tag, aux_appended_type,
                 length_aux_data + 1, aux_data);
  if (sam_write1(file_output, gsDataHdr(gs), new_rec) < 0) {
    fprintf(stderr, "Error: failed to write sam record. \n");
    printSamRecord_brief(gs, new_rec);
    printf("aux data: %s\n", aux_data);
    printf("length aux data: %d\n", length_aux_data);
    exit(EXIT_FAILURE);
  }

  bam_destroy1(new_rec);
  destroy_AlignResult(ar_lpart);
  destroy_AlignResult(ar_rpart);
  return;
}

/**
 * @brief  Find all variants that should be integrated into reference sequence
 * * [lbound, rbound]. Pack them in the form of an ervArray.
 */
static inline void integration_generate_ervArray(
    int64_t lbound, int64_t rbound, const char *rname, GenomeVcf_bplus *gv,
    Element_RecVcf ***ret_ervArray, int *ret_cnt_integrated_variants) {
  assert(lbound <= rbound);
  // printf("generated ervArray lbound: %" PRId64 ", rbound: %" PRId64 "\n",
  //        lbound, rbound);
  // ------- get all variants' combination within the interval -----
  // 1st loop - calculate number of vcf records that needs integration
  int cnt_integrated_variants = 0;
  RecVcf_bplus *rv_tmp = genomeVcf_bplus_getRecAfterPos(gv, rname, lbound);
  RecVcf_bplus *first_rv = rv_tmp;
  while (rv_tmp != NULL) {
    if (count_integrated_allele(rv_tmp, lbound, rbound) > 0) {
      cnt_integrated_variants++;
    } else {
      if (rv_pos(rv_tmp) > rbound) {
        break;
      }
    }
    rv_tmp = next_RecVcf_bplus(rv_tmp);
  }
  // 2nd loop - save pointers to vcf records that needs integration and
  // calculate number of alleles of each vcf record that needs integration.
  Element_RecVcf **ervArray = (Element_RecVcf **)calloc(
      cnt_integrated_variants, sizeof(Element_RecVcf *));
  rv_tmp = first_rv;
  int idx_integrated_variants = 0;
  while (rv_tmp != NULL) {
    int cnt_integrated_allele = count_integrated_allele(rv_tmp, lbound, rbound);
    if (cnt_integrated_allele > 0) {
      // Find all alleles that needs integration on this vcf record
      ervArray[idx_integrated_variants] =
          (Element_RecVcf *)malloc(sizeof(Element_RecVcf));
      ervArray[idx_integrated_variants]->rv = rv_tmp;
      ervArray[idx_integrated_variants]->alleleIdx =
          (int *)calloc(cnt_integrated_allele, sizeof(int));
      ervArray[idx_integrated_variants]->alleleCnt = 0;

      int tmp_idx_allele = 0;
      for (int i = 0; i < rv_alleleCnt(rv_tmp); i++) {
        if (ifCanIntegrateAllele(rv_tmp, i, lbound, rbound) == 1) {
          ervArray[idx_integrated_variants]->alleleIdx[tmp_idx_allele] = i;
          ervArray[idx_integrated_variants]->alleleCnt++;
          tmp_idx_allele++;
        }
      }
      idx_integrated_variants++;
    } else {
      if (rv_pos(rv_tmp) > rbound) {
        break;
      }
    }
    rv_tmp = next_RecVcf_bplus(rv_tmp);
  }

  // printf("rvArray: \n");
  // for (int i = 0; i < cnt_integrated_variants; i++) {
  //   genomeVcf_bplus_printRec(gv, ervArray[i]->rv);
  //   printf("\tintegrated alleles' idxes: ");
  //   for (int j = 0; j < ervArray[i]->alleleCnt; j++) {
  //     printf("%d ", ervArray[i]->alleleIdx[j]);
  //   }
  //   printf("\n");
  // }
  *ret_ervArray = ervArray;
  *ret_cnt_integrated_variants = cnt_integrated_variants;
}

static inline void integration_select_and_integrate(
    Element_RecVcf *ervArray_lpart[], int length_ervArray_lpart,
    Element_RecVcf *ervArray_rpart[], int length_ervArray_rpart,
    int64_t lbound_var, int64_t rbound_var, int64_t lbound_M, int64_t rbound_M,
    RecSam *rec_rs, int64_t id_rec, GenomeFa *gf, GenomeSam *gs,
    GenomeVcf_bplus *gv, samFile *file_output) {
  if (length_ervArray_lpart == 0) {
    // ------------------------ Process right part -----------------------
    int *ervIdxes_rpart = (int *)calloc(length_ervArray_rpart, sizeof(int));
    for (int l = 0; l < length_ervArray_rpart; l++) ervIdxes_rpart[l] = l;
    for (int l = 1; l < length_ervArray_rpart + 1; l++) {
      Combinations *cbs_rpart =
          calculate_combinations(ervIdxes_rpart, length_ervArray_rpart, l);
      // printf("erv combi rpart (size: %d) cnt: %d\n", l, cbs_rpart->cnt);
      for (int m = 0; m < cbs_rpart->cnt; m++) {
        // printf("erv combi rpart [%d]: ", m);
        // for (int n = 0; n < cbs_rpart->cnt; n++) {
        //   printf("%d ", cbs_rpart->combis[m][n]);
        // }
        // printf("\n");
        Combinations_alleles *acbs_rpart = calculate_combinations_alleles(
            ervArray_rpart, length_ervArray_rpart, cbs_rpart->combis[m],
            cbs_rpart->length);
        // printf("allele combi cnt rpart: %d\n", acbs_rpart->cnt);
        // print_combinations_alleles(acbs_rpart);
        for (int n = 0; n < acbs_rpart->cnt; n++) {
          // ------------- Do realignment with left part unmodified ------------
          integration_integrate(
              NULL, NULL, NULL, 0, 0, ervArray_rpart, acbs_rpart->combi_rv,
              acbs_rpart->combis_allele[n], acbs_rpart->length,
              length_ervArray_rpart, lbound_var, rbound_var, lbound_M, rbound_M,
              rec_rs, id_rec, gf, gs, gv, file_output);
        }
        for (int n = 0; n < acbs_rpart->cnt; n++) {
          free(acbs_rpart->combis_allele[n]);
        }
        free(acbs_rpart->combis_allele);
        destroy_combinations_alleles(acbs_rpart);
      }
      for (int m = 0; m < cbs_rpart->cnt; m++) {
        free(cbs_rpart->combis[m]);
      }
      free(cbs_rpart->combis);
      destroy_combinations(cbs_rpart);
    }
    free(ervIdxes_rpart);
    return;
  } else {
    // --------------------------- Process left part ---------------------------
    int *ervIdxes_lpart = (int *)calloc(length_ervArray_lpart, sizeof(int));
    for (int i = 0; i < length_ervArray_lpart; i++) ervIdxes_lpart[i] = i;
    for (int i = 1; i < length_ervArray_lpart + 1; i++) {
      Combinations *cbs_lpart =
          calculate_combinations(ervIdxes_lpart, length_ervArray_lpart, i);
      // printf("erv combi lpart (size: %d) cnt: %d\n", i, cbs_lpart->cnt);
      for (int j = 0; j < cbs_lpart->cnt; j++) {
        // printf("erv combi lpart [%d]: ", j);
        // for (int k = 0; k < cbs_lpart->length; k++) {
        //   printf("%d ", cbs_lpart->combis[j][k]);
        // }
        // printf("\n");
        Combinations_alleles *acbs_lpart = calculate_combinations_alleles(
            ervArray_lpart, length_ervArray_lpart, cbs_lpart->combis[j],
            cbs_lpart->length);
        // printf("allele combi cnt lpart: %d\n", acbs_lpart->cnt);
        // print_combinations_alleles(acbs_lpart);

        for (int k = 0; k < acbs_lpart->cnt; k++) {
          if (length_ervArray_rpart == 0) {
            // ----------- Do realignment with right part unmodified -----------
            integration_integrate(ervArray_lpart, acbs_lpart->combi_rv,
                                  acbs_lpart->combis_allele[k],
                                  acbs_lpart->length, length_ervArray_lpart,
                                  NULL, NULL, NULL, 0, 0, lbound_var,
                                  rbound_var, lbound_M, rbound_M, rec_rs,
                                  id_rec, gf, gs, gv, file_output);
          } else {
            // ----------------------- Process right part ----------------------
            int *ervIdxes_rpart =
                (int *)calloc(length_ervArray_rpart, sizeof(int));
            for (int l = 0; l < length_ervArray_rpart; l++)
              ervIdxes_rpart[l] = l;
            for (int l = 1; l < length_ervArray_rpart + 1; l++) {
              Combinations *cbs_rpart = calculate_combinations(
                  ervIdxes_rpart, length_ervArray_rpart, l);
              // printf("erv combi rpart (size: %d) cnt: %d\n", l,
              // cbs_rpart->cnt);
              for (int m = 0; m < cbs_rpart->cnt; m++) {
                // printf("erv combi rpart [%d]: ", m);
                // for (int n = 0; n < cbs_rpart->cnt; n++) {
                //   printf("%d ", cbs_rpart->combis[m][n]);
                // }
                // printf("\n");
                Combinations_alleles *acbs_rpart =
                    calculate_combinations_alleles(
                        ervArray_rpart, length_ervArray_rpart,
                        cbs_rpart->combis[m], cbs_rpart->length);
                // printf("allele combi cnt rpart: %d\n", acbs_rpart->cnt);
                // print_combinations_alleles(acbs_rpart);
                for (int n = 0; n < acbs_rpart->cnt; n++) {
                  // Do realignment
                  integration_integrate(
                      ervArray_lpart, acbs_lpart->combi_rv,
                      acbs_lpart->combis_allele[k], acbs_lpart->length,
                      length_ervArray_lpart, ervArray_rpart,
                      acbs_rpart->combi_rv, acbs_rpart->combis_allele[n],
                      acbs_rpart->length, length_ervArray_rpart, lbound_var,
                      rbound_var, lbound_M, rbound_M, rec_rs, id_rec, gf, gs,
                      gv, file_output);
                }
                for (int n = 0; n < acbs_rpart->cnt; n++) {
                  free(acbs_rpart->combis_allele[n]);
                }
                free(acbs_rpart->combis_allele);
                destroy_combinations_alleles(acbs_rpart);
              }
              for (int m = 0; m < cbs_rpart->cnt; m++) {
                free(cbs_rpart->combis[m]);
              }
              free(cbs_rpart->combis);
              destroy_combinations(cbs_rpart);
            }
            free(ervIdxes_rpart);
          }
        }
        for (int k = 0; k < acbs_lpart->cnt; k++) {
          free(acbs_lpart->combis_allele[k]);
        }
        free(acbs_lpart->combis_allele);
        destroy_combinations_alleles(acbs_lpart);
      }
      for (int j = 0; j < cbs_lpart->cnt; j++) {
        free(cbs_lpart->combis[j]);
      }
      free(cbs_lpart->combis);
      destroy_combinations(cbs_lpart);
    }
    free(ervIdxes_lpart);
    return;
  }
}

typedef struct _define_ThreadArgs {
  int64_t id;  // identifier for the thread
  Options *opts;
  GenomeFa *gf;
  GenomeSam *gs;
  GenomeVcf_bplus *gv;
  int64_t id_sam_start;  // index for sam record when the thread starts
  int64_t id_sam_end;    // index for sam record when the thread ends
} ThreadArgs;

void *integration_threads(void *args) {
  ThreadArgs *args_thread = (ThreadArgs *)args;

  // printf("thread (%" PRId64 ") process sam rec from %" PRId64 " to %" PRId64
  //        "\n",
  //        args_thread->id, args_thread->id_sam_start,
  //        args_thread->id_sam_end);
  // Execute integrations process

  char path_outputFile[strlen(getOutputFile(args_thread->opts)) + 64];
  sprintf(path_outputFile, "%s.thread%" PRId64 "",
          getOutputFile(args_thread->opts), args_thread->id);

  // printf("output file name for thread (%" PRId64 "): %s\n", args_thread->id,
  //        path_outputFile);

  samFile *file_output = NULL;
  file_output = sam_open(path_outputFile, "w");
  if (file_output == NULL) {
    fprintf(stderr, "Error: cannot open file %s with mode \"w\"\n",
            path_outputFile);
    exit(EXIT_FAILURE);
  }

  // Write header into file
  sam_hdr_t *header = gsDataHdr(args_thread->gs);
  if (sam_hdr_write(file_output, header) < 0) {
    fprintf(stderr, "Error: failed writing header for %s\n", path_outputFile);
    exit(EXIT_FAILURE);
  }

  // Extract arguments from thread input
  GenomeFa *gf = args_thread->gf;
  GenomeSam *gs = args_thread->gs;
  GenomeVcf_bplus *gv = args_thread->gv;

  // Iterate all sam records and locate their corresponding variants.
  GenomeSamIterator *gsIt = init_GenomeSamIterator(gs);
  ChromSam *cs_tmp = gsItNextChrom(gsIt);
  RecSam *rs_tmp = gsItNextRec(gsIt);

  int64_t id_rec = 0;  // Id of record processed in this thread
  int64_t id_rec_start = args_thread->id_sam_start;  // included
  int64_t id_rec_end = args_thread->id_sam_end;      // included
  while (rs_tmp != NULL) {
    // -------- only handle records within [id_start, id_end] --------
    if (id_rec < id_rec_start) {
      rs_tmp = gsItNextRec(gsIt);
      if (rs_tmp == NULL) {
        cs_tmp = gsItNextChrom(gsIt);
        rs_tmp = gsItNextRec(gsIt);
      }
      id_rec++;
      continue;
    }
    if (id_rec > id_rec_end) {
      break;
    }
    // ---------- get information of temporary sam record ------------
    const char *rname_read = rsDataRname(gs, rs_tmp);
    int64_t lbound_read = rsDataPos(rs_tmp);  // 1-based, included
    int64_t rbound_read = lbound_read - 1;    // 1-based, included

    // ------------- find the longest 'M' area in cigar --------------
    int64_t tmp_lbound = 1;
    int64_t lbound_M_ref = 1;  // 1-based, included
    uint32_t cnt_cigar = rs_cigar_cnt(rs_tmp);
    uint32_t length_M_area = 0;
    for (int i = 0; i < cnt_cigar; i++) {
      uint32_t tmp_length = rs_cigar_oplen(rs_tmp, i);
      char cigar_opChar = rs_cigar_opChar(rs_tmp, i);
      if (cigar_opChar == 'M') {
        if (tmp_length >= length_M_area) {
          length_M_area = tmp_length;
          lbound_M_ref = tmp_lbound;
        }
      }
      if (cigar_opChar == 'D') {
        rbound_read += tmp_length;
        tmp_lbound += tmp_length;
      } else if (cigar_opChar == 'I') {
        // rbound_read will not move on the reference in such case
      } else if (cigar_opChar == 'P' || cigar_opChar == 'N') {
        length_M_area = 0;
        break;
        // Ignore records containing 'P' and 'N'
      } else {  // For cigarOp 'MX=SH', move lbound_M_ref on the reference
        rbound_read += tmp_length;
        tmp_lbound += tmp_length;
      }
    }
    int64_t rbound_M_ref =
        lbound_M_ref + length_M_area - 1;  // 1-based, included
    // Add offset based on the start position of read
    lbound_M_ref += lbound_read - 1;
    rbound_M_ref += lbound_read - 1;
    // printf("lbound read: %" PRId64 ", rbound read: %" PRId64 "\n",
    // lbound_read,
    //        rbound_read);
    // printf("length M area: %" PRIu32 ", lbound cigar M: % " PRId64
    //        ", rbound cigar M: %" PRId64 "\n",
    //        length_M_area, lbound_M_ref, rbound_M_ref);

    // ----------- ignore reads as following in this program ---------
    // empty rname, empty cigar, or no M_area
    if (rname_read == NULL || length_M_area == 0) {
      rs_tmp = gsItNextRec(gsIt);
      if (rs_tmp == NULL) {
        cs_tmp = gsItNextChrom(gsIt);
        rs_tmp = gsItNextRec(gsIt);
      }
      continue;
    }

    // -------- get boundaries of area for selecting variants --------
    int64_t lbound_variant = 0;            // 1-based, included
    int64_t rbound_variant = rbound_read;  // 1-based, included
    switch (opt_integration_strategy(args_thread->opts)) {
      case _OPT_INTEGRATION_SNPONLY: {
        lbound_variant = lbound_read - integration_sv_min_len;
        break;
      }
      case _OPT_INTEGRATION_SVONLY:
      case _OPT_INTEGRATION_ALL: {
        lbound_variant = lbound_read - integration_sv_max_len;
        break;
      }
      default: {
        fprintf(stderr, "Error: no such strategy for integration.\n");
        exit(EXIT_FAILURE);
      }
    }
    if (lbound_variant <= 0) lbound_variant = 1;
    // printf("lbound variant: %" PRId64 ", rbound variant: %" PRId64 "\n",
    //        lbound_variant, rbound_variant);

    // Split the area into 2 parts
    // ** lbound_var **1** lbound_M_ref M..M rbound_M_ref **2*** rbound_var **

    Element_RecVcf **ervArray_lpart = NULL;
    Element_RecVcf **ervArray_rpart = NULL;
    int cnt_integrated_variants_lpart = 0;
    int cnt_integrated_variants_rpart = 0;
    // printf("L part generated ervArray - ");
    integration_generate_ervArray(lbound_variant, lbound_M_ref, rname_read, gv,
                                  &ervArray_lpart,
                                  &cnt_integrated_variants_lpart);
    // printf("R part generated ervArray - ");
    integration_generate_ervArray(rbound_M_ref, rbound_variant, rname_read, gv,
                                  &ervArray_rpart,
                                  &cnt_integrated_variants_rpart);

    // ---------------- select alleles and integrate -----------------
    integration_select_and_integrate(
        ervArray_lpart, cnt_integrated_variants_lpart, ervArray_rpart,
        cnt_integrated_variants_rpart, lbound_variant, rbound_variant,
        lbound_M_ref, rbound_M_ref, rs_tmp, id_rec, gf, gs, gv, file_output);

    // ----------------------- free memories -------------------------
    for (int i = 0; i < cnt_integrated_variants_lpart; i++) {
      free(ervArray_lpart[i]->alleleIdx);
      free(ervArray_lpart[i]);
    }
    free(ervArray_lpart);
    for (int i = 0; i < cnt_integrated_variants_rpart; i++) {
      free(ervArray_rpart[i]->alleleIdx);
      free(ervArray_rpart[i]);
    }
    free(ervArray_rpart);
    // printSamRecord_brief(gs, rsData(rs_tmp));
    // printf("*****************************************************\n");

    // --------------------- keep on iterating -----------------------
    rs_tmp = gsItNextRec(gsIt);
    if (rs_tmp == NULL) {
      cs_tmp = gsItNextChrom(gsIt);
      rs_tmp = gsItNextRec(gsIt);
    }
    id_rec++;
  }

  destroy_GenomeSamIterator(gsIt);

  sam_close(file_output);

  return (void *)(args_thread->id);
}

void integration(Options *opts) {
  check_files(opts);

  // Init alignment parameters
  alignInitialize(getMatch(opts), getMismatch(opts), getGapopen(opts),
                  getGapextension(opts));
  // Init paramters for ksw2 specially
  alignInitialize_ksw2(KSW2_DEFAULT_BANDWIDTH, KSW2_DEFAULT_ZDROP,
                       KSW2_DEFAULT_FLAG);
  integration_strategy = opt_integration_strategy(opts);
  integration_sv_min_len = getSVminLen(opts);
  integration_sv_max_len = getSVmaxLen(opts);

  clock_t time_start = 0;
  clock_t time_end = 0;
  // Init structures (data storage and access)
  time_start = clock();
  GenomeFa *gf = init_GenomeFa();
  loadGenomeFaFromFile(gf, getFaFile(opts));
  time_end = clock();
  printf("... %s loaded. time: %fs\n", getFaFile(opts),
         time_convert_clock2second(time_start, time_end));
  time_start = clock();
  GenomeSam *gs = init_GenomeSam();
  loadGenomeSamFromFile(gs, getSamFile(opts));
  time_end = clock();
  printf("... %s loaded. time: %fs\n", getSamFile(opts),
         time_convert_clock2second(time_start, time_end));
  time_start = clock();
  GenomeVcf_bplus *gv = genomeVcf_bplus_loadFile(getVcfFile(opts), 7, 6);
  time_end = clock();
  printf("... %s loaded. time: %fs\n", getVcfFile(opts),
         time_convert_clock2second(time_start, time_end));

  time_start = clock();
  const int cnt_thread = opt_threads(opts) > 0 ? opt_threads(opts) : 1;
  const int64_t cnt_rec_sam = gsDataRecCnt(gs);
  pthread_t threads[cnt_thread];
  ThreadArgs args_thread[cnt_thread];

  // Assign arguments for threads
  for (int i = 0; i < cnt_thread; i++) {
    args_thread[i].id = i;
    args_thread[i].opts = opts;
    args_thread[i].gf = gf;
    args_thread[i].gs = gs;
    args_thread[i].gv = gv;
    args_thread[i].id_sam_start = (cnt_rec_sam / cnt_thread) * i;
  }
  for (int i = 0; i < cnt_thread - 1; i++) {
    args_thread[i].id_sam_end = args_thread[i + 1].id_sam_start - 1;
  }
  args_thread[cnt_thread - 1].id_sam_end = cnt_rec_sam - 1;

  // Create threads
  for (int i = 0; i < cnt_thread; i++) {
    if (pthread_create(&threads[i], NULL, integration_threads,
                       (void *)&args_thread[i]) != 0) {
      fprintf(
          stderr,
          "Error: failed to create thread to handle sam records from %" PRId64
          " to %" PRId64 "\n",
          args_thread->id_sam_start, args_thread->id_sam_end);
      exit(EXIT_FAILURE);
    }
  }

  // End threads
  for (int i = 0; i < cnt_thread; i++) {
    void *thread_ret = 0;
    pthread_join(threads[i], &thread_ret);
    printf("thread (%" PRId64 ") ended\n", (int64_t)thread_ret);
  }
  time_end = clock();
  printf("... integration finished. Total time: %fs\n",
         time_convert_clock2second(time_start, time_end));

  // Free structures

  destroy_GenomeFa(gf);
  destroy_GenomeSam(gs);
  destroy_GenomeVcf_bplus(gv);
  return;
}
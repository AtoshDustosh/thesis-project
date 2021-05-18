#include "integrateVcfToSam.h"

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
  // This is based on the assumption that varEndPos > varStartPoos
  int64_t varStartPos = rv_pos(rv);
  int64_t varAffectedLength = rv_alleleCoverLength(rv, alleleIdx);
  int64_t varEndPos = varStartPos + varAffectedLength - 1;
  return varStartPos <= endPos && varEndPos >= startPos;
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
        assert((var_alt[0] != '<') ||
               (fprintf(stderr,
                        "Error: variants with tags should be ignored when "
                        "loading.\n") < 0));
        int length_ref = strlen(var_ref);
        int length_alt = strlen(var_alt);
        if (length_ref == 1) {  // Insertion
          // Insertion doesn't affect the following variants
          if (varStartPos >= startPos && varStartPos <= endPos)
            cnt_integrated_allele++;
        } else if (length_alt == 1) {  // Deletion
          // Deletion might affect the following variants
          int64_t varEndPos = varStartPos + length_ref - 1;
          if (varEndPos >= startPos && varStartPos <= endPos)
            cnt_integrated_allele++;
        } else {
          // For those records like (REF ALT) = (ACG A,ACGTT) or like (REF ALT)
          // = (AG CT)
          int64_t varEndPos = varStartPos + length_ref - 1;
          if (varEndPos >= startPos && varStartPos <= endPos)
            cnt_integrated_allele++;
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

static inline void integration_integrate(
    Element_RecVcf *ervArray_lpart[], int ervCombi_lpart[],
    int alleleCombi_lpart[], int length_combi_lpart,
    Element_RecVcf *ervArray_rpart[], int ervCombi_rpart[],
    int alleleCombi_rpart[], int length_combi_rpart, int64_t lbound_var,
    int64_t rbound_var, int64_t lbound_M, int64_t rbound_M, RecSam *rec_rs,
    GenomeFa *gf, GenomeSam *gs, GenomeVcf_bplus *gv, samFile *file_output) {
  // Print information collected for this integration task
  printf("lbound var: %" PRId64 ", rbound var: %" PRId64 "\n", lbound_var,
         rbound_var);
  printf("lbound_M: %" PRId64 ", rbound_M: %" PRId64 "\n", lbound_M, rbound_M);
  printf("length combi lpart: %d, combi rpart: %d\n", length_combi_lpart,
         length_combi_rpart);
  // printf("ervArray (L): \n");
  // for (int i = 0; i < length_combi_lpart; i++) {
  //   genomeVcf_bplus_printRec(gv, ervArray_lpart[i]->rv);
  // }
  // printf("\n");
  // printf("(rv, allele) selected (L): \n");
  // for (int i = 0; i < length_combi_lpart; i++) {
  //   printf("(%d,%d) ", ervCombi_lpart[i], alleleCombi_lpart[i]);
  // }
  // printf("\n");
  // printf("ervArray (R): \n");
  // for (int i = 0; i < length_combi_rpart; i++) {
  //   genomeVcf_bplus_printRec(gv, ervArray_rpart[i]->rv);
  // }
  // printf("\n");
  // printf("(rv, allele) selected (R): \n");
  // for (int i = 0; i < length_combi_rpart; i++) {
  //   printf("(%d,%d) ", ervCombi_rpart[i], alleleCombi_rpart[i]);
  // }
  // printf("\n");
  // printf("*****************************************************\n");

  // Integrate the right part
  // TODO
  // Integrate the left part
  // TODO
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
    RecSam *rec_rs, GenomeFa *gf, GenomeSam *gs, GenomeVcf_bplus *gv,
    samFile *file_output) {
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
              NULL, NULL, NULL, 0, ervArray_rpart, acbs_rpart->combi_rv,
              acbs_rpart->combis_allele[n], acbs_rpart->length, lbound_var,
              rbound_var, lbound_M, rbound_M, rec_rs, gf, gs, gv, file_output);
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
                                  acbs_lpart->length, NULL, NULL, NULL, 0,
                                  lbound_var, rbound_var, lbound_M, rbound_M,
                                  rec_rs, gf, gs, gv, file_output);
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
                      ervArray_rpart, acbs_rpart->combi_rv,
                      acbs_rpart->combis_allele[n], acbs_rpart->length,
                      lbound_var, rbound_var, lbound_M, rbound_M, rec_rs, gf,
                      gs, gv, file_output);
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

  printf("output file name for thread (%" PRId64 "): %s\n", args_thread->id,
         path_outputFile);

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

  const int sv_min_len = getSVminLen(args_thread->opts);
  const int sv_max_len = getSVmaxLen(args_thread->opts);

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
    printSamRecord_brief(gs, rsData(rs_tmp));
    const char *rname_read = rsDataRname(gs, rs_tmp);
    int64_t lbound_read = rsDataPos(rs_tmp);  // 1-based, included
    uint32_t length_read = rsDataSeqLength(rs_tmp);
    int64_t rbound_read = lbound_read;  // 1-based, included

    // ------------- find the longest 'M' area in cigar --------------
    int64_t tmp_lbound = 0;
    int64_t lbound_M = 0;  // 1-based, included
    uint32_t cnt_cigar = rs_cigar_cnt(rs_tmp);
    uint32_t length_M_area = 0;
    for (int i = 0; i < cnt_cigar; i++) {
      uint32_t tmp_length = rs_cigar_oplen(rs_tmp, i);
      char cigar_opChar = rs_cigar_opChar(rs_tmp, i);
      if (cigar_opChar == 'M') {
        if (tmp_length >= length_M_area) {
          length_M_area = tmp_length;
          lbound_M = tmp_lbound;
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
      } else {  // For cigarOp 'MX=SH', move lbound_M on the reference
        rbound_read += tmp_length;
        tmp_lbound += tmp_length;
      }
    }
    int64_t rbound_M = lbound_M + length_M_area - 1;  // 1-based, included
    // Add offset based on the start position of read
    lbound_M += lbound_read - 1;
    rbound_M += lbound_read - 1;
    // printf("lbound read: %" PRId64 ", rbound read: %" PRId64 "\n",
    // lbound_read,
    //        rbound_read);
    // printf("length M area: %" PRIu32 ", lbound cigar M: % " PRId64
    //        ", rbound cigar M: %" PRId64 "\n",
    //        length_M_area, lbound_M, rbound_M);

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
        lbound_variant = lbound_read - sv_min_len;
        break;
      }
      case _OPT_INTEGRATION_SVONLY:
      case _OPT_INTEGRATION_ALL: {
        lbound_variant = lbound_read - sv_max_len;
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
    // ** lbound_var **1** lbound_M M..M rbound_M **2*** rbound_var **

    Element_RecVcf **ervArray_lpart = NULL;
    Element_RecVcf **ervArray_rpart = NULL;
    int cnt_integrated_variants_lpart = 0;
    int cnt_integrated_variants_rpart = 0;
    // printf("L part generated ervArray - ");
    integration_generate_ervArray(lbound_variant, lbound_M, rname_read, gv,
                                  &ervArray_lpart,
                                  &cnt_integrated_variants_lpart);
    // printf("R part generated ervArray - ");
    integration_generate_ervArray(rbound_M, rbound_variant, rname_read, gv,
                                  &ervArray_rpart,
                                  &cnt_integrated_variants_rpart);

    // ---------------- select alleles and integrate -----------------
    integration_select_and_integrate(
        ervArray_lpart, cnt_integrated_variants_lpart, ervArray_rpart,
        cnt_integrated_variants_rpart, lbound_variant, rbound_variant, lbound_M,
        rbound_M, rs_tmp, gf, gs, gv, file_output);

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

  // Init structures (data storage and access)
  GenomeFa *gf = init_GenomeFa();
  loadGenomeFaFromFile(gf, getFaFile(opts));
  GenomeSam *gs = init_GenomeSam();
  loadGenomeSamFromFile(gs, getSamFile(opts));
  GenomeVcf_bplus *gv = genomeVcf_bplus_loadFile(getVcfFile(opts), 7, 6);

  const int cnt_thread = opt_threads(opts);
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

  // Free structures

  destroy_GenomeFa(gf);
  destroy_GenomeSam(gs);
  destroy_GenomeVcf_bplus(gv);
  return;
}
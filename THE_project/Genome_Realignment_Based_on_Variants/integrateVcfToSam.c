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

static inline int integrate_realign_output(Element_RecVcf *ervArray[],
                                           int *ervCombi, int *combi_allele,
                                           int length_combi, GenomeFa *gf,
                                           GenomeSam *gs, GenomeVcf_bplus *gv,
                                           RecSam *rs, samFile *file_output) {
  // TODO Check selected variants and extract ref sequence. If there exists DEL,
  // you might need to extract longer ref sequence.
  return 0;
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
    const char *qname_read = rsDataQname(rs_tmp);
    int64_t lbound_read = rsDataPos(rs_tmp);  // included
    uint32_t length_read = rsDataSeqLength(rs_tmp);
    char *seq_read = rsDataSeq(rs_tmp);

    int64_t rbound_read = lbound_read + length_read - 1;

    // ---------- ignore those unmapped reads in this program --------
    if (rname_read == NULL) {
      rs_tmp = gsItNextRec(gsIt);
      if (rs_tmp == NULL) {
        cs_tmp = gsItNextChrom(gsIt);
        rs_tmp = gsItNextRec(gsIt);
      }
      continue;
    }

    // -------- get boundaries of area for selecting variants --------
    int64_t lbound_variant = 0;
    int64_t rbound_variant = rbound_read;
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
    // printf("lbound: %" PRId64 ", rbound: %" PRId64 "\n", lbound_variant,
    //        rbound_variant);

    // ------- get all variants' combination within the interval -----
    // 1st loop - calculate number of vcf records that needs integration
    int cnt_integrated_variants = 0;
    RecVcf_bplus *rv_tmp =
        genomeVcf_bplus_getRecAfterPos(gv, rname_read, lbound_variant);
    RecVcf_bplus *first_rv = rv_tmp;
    while (rv_tmp != NULL) {
      if (count_integrated_allele(rv_tmp, lbound_variant, rbound_variant) > 0) {
        cnt_integrated_variants++;
      } else {
        if (rv_pos(rv_tmp) > rbound_variant) {
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
      int cnt_integrated_allele =
          count_integrated_allele(rv_tmp, lbound_variant, rbound_variant);
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
          if (ifCanIntegrateAllele(rv_tmp, i, lbound_variant, rbound_variant) ==
              1) {
            ervArray[idx_integrated_variants]->alleleIdx[tmp_idx_allele] = i;
            ervArray[idx_integrated_variants]->alleleCnt++;
            tmp_idx_allele++;
          }
        }
        idx_integrated_variants++;
      } else {
        if (rv_pos(rv_tmp) > rbound_variant) {
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

    // -------------------- perform integration ----------------------
    int *ervIdxes = (int *)calloc(cnt_integrated_variants, sizeof(int));
    for (int i = 0; i < cnt_integrated_variants; i++) ervIdxes[i] = i;
    int cnt_new_recsam = 0;
    for (int i = 1; i < cnt_integrated_variants + 1; i++) {
      Combinations *cbs =
          calculate_combinations(ervIdxes, cnt_integrated_variants, i);
      // printf("erv combi(size: %d) cnt: %d\n", i, cbs->cnt);
      for (int j = 0; j < cbs->cnt; j++) {
        // printf("erv combi[%d]: ", j);
        // for (int k = 0; k < cbs->length; k++) {
        //   printf("%d ", cbs->combis[j][k]);
        // }
        // printf("\n");
        Combinations_alleles *acbs = calculate_combinations_alleles(
            ervArray, cnt_integrated_variants, cbs->combis[j], cbs->length);
        // printf("allele combi cnt: %d\n", acbs->cnt);
        // print_combinations_alleles(acbs);
        for (int k = 0; k < acbs->cnt; k++) {
          cnt_new_recsam += integrate_realign_output(
              ervArray, acbs->combi_rv, acbs->combis_allele[k], acbs->length,
              gf, gs, gv, rs_tmp, file_output);
        }
        destroy_combinations_alleles(acbs);
      }
      destroy_combinations(cbs);
    }

    // ----------------------- free memories -------------------------
    free(seq_read);
    free(ervIdxes);
    for (int i = 0; i < cnt_integrated_variants; i++) {
      free(ervArray[i]->alleleIdx);
      free(ervArray[i]);
    }
    free(ervArray);
    printf("*****************************************************\n");

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
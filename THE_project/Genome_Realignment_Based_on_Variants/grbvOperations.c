#include "grbvOperations.h"

/**
 * @brief  An auxiliary data structure designed for pointint to alleles among an
 * array of RecVcf from which the combinations of alleles are selected. Used
 * with RecVcf[].
 */
typedef struct _define_PendingAlleles {
  int rvIdx;
  /**
   * @brief  Although most vcf records only have one allele, just in case of the
   * special ones, we use an array for the alleles' indexes in a vcf record.
   */
  int *alleleIdx;
  int alleleCnt;
} PendingAlleles;

//  --------------------------- split line ---------------------------
//  --------------------------- split line ---------------------------
//  --------------------------- split line ---------------------------
//  --------------------------- split line ---------------------------
//  --------------------------- split line ---------------------------

/**
 * @brief NOT RECOMMENDED. Use method "countIntegratedAllele > 0" instead. Judge
 * whether a variant could be integrated or not.
 *
 * @param rv variant
 * @param startPos 1-based position of the start point of refSeq
 * @param endPos 1-based position of the end point of refSeq
 * @return int 1 if the variant could be integrated; 0 otherwise.
 */
static inline int ifCanIntegrateVar(RecVcf *rv, int64_t startPos,
                                    int64_t endPos) {
  // this is based on the assumption that varEndPos > varStartPos
  int64_t varStartPos = rvDataPos(rv);
  // calculate the affected length
  int64_t varAffectedLength = 0;
  for (int i = 0; i < rvDataAlleleCnt(rv); i++) {
    switch (bcf_get_variant_type(rvData(rv), i)) {
      case VCF_REF: {
        break;
      }
      case VCF_SNP: {
        break;
      }
      case VCF_INDEL: {
        const char *varRef = rvData(rv)->d.allele[0];
        const char *varAlt = rvData(rv)->d.allele[i];
        int refLength = strlen(varRef);
        int altLength = strlen(varAlt);
        if (refLength == 1) {  // INS
          // insertion doesn't affect the following varaints
        } else if (altLength == 1) {  // DEL
          // deletioin may affect the following variants
          if (varAffectedLength < refLength) varAffectedLength = refLength - 1;
        } else {
          // for those records like this, (REF, ALT) = (ACG, A,ACGTT)
          // or like this, (REF, ALT) = (AG, CT)
          if (varAffectedLength < refLength) varAffectedLength = refLength - 1;
        }
        break;
      }
      default: {
        // ignore other kinds of varaints
      }
    }
  }
  int64_t varEndPos = varStartPos + varAffectedLength;
  return varStartPos <= endPos && varEndPos >= startPos;
}

/**
 * @brief  This method is used to count the alleles that should be integrated
 * within a region.
 */
static inline int countIntegratedAllele(RecVcf *rv, int64_t startPos,
                                        int64_t endPos) {
  // this is based on the assumption that varEndPos > varStartPos
  int64_t varStartPos = rvDataPos(rv);
  int integratedAlleleCnt = 0;
  // calculate the affected length
  for (int i = 0; i < rvDataAlleleCnt(rv); i++) {
    switch (bcf_get_variant_type(rvData(rv), i)) {
      case VCF_REF: {
        break;
      }
      case VCF_SNP: {
        if (varStartPos >= startPos && varStartPos <= endPos)
          integratedAlleleCnt++;
        break;
      }
      case VCF_INDEL: {
        const char *varRef = rvData(rv)->d.allele[0];
        const char *varAlt = rvData(rv)->d.allele[i];
        int refLength = strlen(varRef);
        int altLength = strlen(varAlt);
        if (refLength == 1) {  // INS
          // insertion doesn't affect the following varaints
          if (varStartPos >= startPos && varStartPos <= endPos)
            integratedAlleleCnt++;
        } else if (altLength == 1) {  // DEL
          // deletioin may affect the following variants
          int64_t varEndPos = varStartPos + refLength - 1;
          if (varEndPos >= startPos && varStartPos <= endPos)
            integratedAlleleCnt++;
        } else {
          // for those records like this, (REF, ALT) = (ACG, A,ACGTT)
          // or like this, (REF, ALT) = (AG, CT)
          int64_t varEndPos = varStartPos + refLength - 1;
          if (varEndPos >= startPos && varStartPos <= endPos)
            integratedAlleleCnt++;
        }
        break;
      }
      default: {
        fprintf(stderr,
                "Warning: variant ignored. varRef: %s, varAlt[%d]: %s\n",
                rvData(rv)->d.allele[0], i, rvData(rv)->d.allele[i]);
        // ignore other kinds of varaints
      }
    }
  }
  return integratedAlleleCnt;
}

static inline int ifCanIntegrateAllele(RecVcf *rv, int alleleIdx, int startPos,
                                       int endPos) {
  // this is based on the assumption that varEndPos > varStartPos
  int64_t varStartPos = rvDataPos(rv);
  // calculate the affected length
  int64_t varAffectedLength = 0;
  switch (bcf_get_variant_type(rvData(rv), alleleIdx)) {
    case VCF_REF: {
      return 0;
      break;
    }
    case VCF_SNP: {
      break;
    }
    case VCF_INDEL: {
      const char *varRef = rvData(rv)->d.allele[0];
      const char *varAlt = rvData(rv)->d.allele[alleleIdx];
      int refLength = strlen(varRef);
      int altLength = strlen(varAlt);
      if (refLength == 1) {  // INS
        // insertion doesn't affect the following varaints
      } else if (altLength == 1) {  // DEL
        // deletioin may affect the following variants
        if (varAffectedLength < refLength) varAffectedLength = refLength - 1;
      } else {
        // for those records like this, (REF, ALT) = (ACG, A,ACGTT)
        // or like this, (REF, ALT) = (AG, CT)
        if (varAffectedLength < refLength) varAffectedLength = refLength - 1;
      }
      break;
    }
    default: {
      // ignore other kinds of varaints
    }
  }
  int64_t varEndPos = varStartPos + varAffectedLength;
  return varStartPos <= endPos && varEndPos >= startPos;
}

/**
 * @brief Find the first vcf record to be integrated.
 * // TODO O(n) time complexity. Bad. Very bad.
 *
 * @param cv ChromVcf object
 * @param startPos 1-based position of the start point of refSeq
 * @param endPos 1-based position of the end point of refSeq
 * @return RecVcf* the first vcf record to be integrated
 */
static inline RecVcf *findFirstVarThatCanBeIntegrated(ChromVcf *cv,
                                                      int64_t startPos,
                                                      int64_t endPos) {
  assert(cv != NULL);
  ChromVcfIterator *cvIt = init_ChromVcfIterator(cv);
  RecVcf *rv = cvItNextRec(cvIt);
  while (rv != NULL) {
    if (ifCanIntegrateVar(rv, startPos, endPos)) {
      destroy_ChromVcfIterator(cvIt);
      return rv;
    }
    rv = cvItNextRec(cvIt);
  }

  destroy_ChromVcfIterator(cvIt);
  return NULL;
}

static inline void integrateVarAndRealign_refactored(
    ElementRecVcf *ervArray[], int *ervCombi, int *alleleCombi, int combiSize,
    int64_t refStartPos, const char *oldRefSeq, const char *readSeq,
    GenomeVcf *gv, RecSam *rs, samFile *fp, GenomeSam *gs) {
  // printf("| integrated alleles: ");
  // for (int i = 0; i < combiSize; i++) {
  //   printf("(%d,%d) ", ervCombi[i], alleleCombi[i]);
  // }
  // printf("\n");
  const int oldRefSeqLen = strlen(oldRefSeq);
  // A buffer for constructing the new reference sequence
  // TODO this buffer cannot support the integration of a large SV
  char buf[oldRefSeqLen * 2];
  memset(buf, 0, oldRefSeqLen * 2);

  // Go along both the oldRefSeq and those alleles.
  int idx_oldRefSeq = 0;
  int idx_newRefSeq = 0;
  // Fix value of start pos, which may be changed by a partially integrated DEL
  // at the start of the region
  int fixPosStart = refStartPos;

  // Start integration
  for (int i = 0; i < combiSize; i++) {
    // Synchronize the positions of buf and allele
    int alleleStartPos = rvDataPos(ervArray[ervCombi[i]]->rv);
    const char *allele_ref = rvDataAllele(ervArray[ervCombi[i]]->rv, 0);
    const char *allele_alt =
        rvDataAllele(ervArray[ervCombi[i]]->rv, alleleCombi[i]);
    int len_allele_ref = strlen(allele_ref);
    int len_allele_alt = strlen(allele_alt);
    int idx_allele_alt = 0;
    // Synchronize allele before integration, especially when the
    // allele is a DEL and the start pos is in front of the region's start
    // pos. Ignore the bases in the front of the allele.
    if (alleleStartPos < refStartPos + idx_oldRefSeq) {
      // fprintf(stderr, "-----------Warning: allele partially integrated.\n");
      idx_allele_alt += refStartPos + idx_oldRefSeq - alleleStartPos;
      len_allele_ref -= refStartPos + idx_oldRefSeq - alleleStartPos;
      fixPosStart += len_allele_ref;
    }
    while (alleleStartPos > refStartPos + idx_oldRefSeq) {
      buf[idx_newRefSeq++] = oldRefSeq[idx_oldRefSeq++];
    }
    // Ignore the bases on old ref according to REF of the allele
    for (int j = 0; j < len_allele_ref; j++) {
      if (idx_oldRefSeq < oldRefSeqLen) {
        idx_oldRefSeq++;
      }
    }
    // Copy the ALT bases of the allele to the new sequence
    while (idx_allele_alt < len_allele_alt) {
      buf[idx_newRefSeq++] = allele_alt[idx_allele_alt++];
    }
    if (idx_oldRefSeq > oldRefSeqLen) {
      fprintf(stderr, "Error: array out of boundary. \n");
      printf("old ref start pos: %" PRId64 ", old ref len: %d\n", refStartPos,
             oldRefSeqLen);
      printf("read seq: %s\n", readSeq);
      printSamRecord_brief(gs, rsData(rs));
      printf("variants: \n");
      for (int i = 0; i < combiSize; i++) {
        printVcfRecord_brief(gv, rvData(ervArray[i]->rv));
      }
      printf("\n");
      printf("var combi: ");
      for (int i = 0; i < combiSize; i++) {
        printf("%d ", ervCombi[i]);
      }
      printf("\n");
      printf("allele combi: ");
      for (int i = 0; i < combiSize; i++) {
        printf("%d ", alleleCombi[i]);
      }
      printf("\n");
      exit(EXIT_FAILURE);
    }
  }
  // Pad the rest unloaded old seq bases
  while (idx_oldRefSeq < oldRefSeqLen)
    buf[idx_newRefSeq++] = oldRefSeq[idx_oldRefSeq++];

  // Calculate new MAPQ
  char *newRefSeq = strdup(buf);
  // printf("oldSeq: %s\n", oldRefSeq);
  // printf("newSeq: %s\n", newRefSeq);
  // printf("readSeq: %s\n", readSeq);
  AlignResult *ar = init_AlignResult();
  align_ssw(newRefSeq, strlen(newRefSeq), readSeq, strlen(readSeq), ar);
  fixPos_AlignResult(ar, fixPosStart);
  // print_AlignResult(ar);

  // Write the result into the output file
  bam1_t *newRec =
      bamSetPosCigarMapq(rsData(rs), arDataPos(ar), arDataReadBegin(ar),
                         arDataReadEnd(ar), arDataCigar(ar), arDataMapQ(ar));
  if (sam_write1(fp, gsDataHdr(gs), newRec) < 0) {
    fprintf(stderr, "Error: failed to write sam record - \n");
    printSamRecord_brief(gs, newRec);
    exit(EXIT_FAILURE);
  }

  free(newRefSeq);
  destroy_AlignResult(ar);
  bam_destroy1(newRec);
  return;
}

void integrateVcfToSam_refactored(Options *opts) {
  if (getSamFile(opts) == NULL || getVcfFile(opts) == NULL ||
      getFaFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments not complete for variants integration - lack "
            "input files *.sam , *.fa or *.vcf. \n");
    exit(EXIT_FAILURE);
  }
  if (getOutputFile(opts) == NULL) {
    setOutputFile(opts, "data/defaultOutput.txt");
  }

  GenomeFa *gf = init_GenomeFa();
  GenomeSam *gs = init_GenomeSam();
  GenomeVcf *gv = init_GenomeVcf();

  loadGenomeFaFromFile(gf, getFaFile(opts));
  loadGenomeSamFromFile(gs, getSamFile(opts));
  loadGenomeVcfFromFile(gv, getVcfFile(opts));

  // Write header for output file (new sam file)
  samFile *op_file = NULL;
  op_file = sam_open(getOutputFile(opts), "w");
  if (op_file == NULL) {
    fprintf(stderr, "Error: cannot open file %s with mode \"w\"\n",
            getOutputFile(opts));
    exit(EXIT_FAILURE);
  }

  sam_hdr_t *header = gsDataHdr(gs);
  // sam_hdr_add_line(header, "SQ", "SN", "ref3", "LN", "5003", NULL);

  if (sam_hdr_write(op_file, header) < 0) {
    fprintf(stderr, "Error: failed to write sam file header. \n");
    exit(EXIT_FAILURE);
  }

  // Iterate all sam records and locate their corresponding variants.
  GenomeSamIterator *gsIt = init_GenomeSamIterator(gs);
  ChromSam *tmpCs = gsItNextChrom(gsIt);
  RecSam *tmpRs = gsItNextRec(gsIt);

  GenomeVcfIterator *gvIt = init_GenomeVcfIterator(gv);
  ChromVcf *tmpCv = NULL;
  RecVcf *tmpRv = NULL;

  const uint8_t ecLen = 5;  // extension length for a read on the ref
  while (tmpRs != NULL) {
    // ------------- get information of temporary sam record ------------
    const char *readRname = rsDataRname(gs, tmpRs);
    const char *readQname = bam_get_qname(rsData(tmpRs));
    int64_t readStartPos = rsDataPos(tmpRs);
    uint32_t readLength = rsDataSeqLength(tmpRs);
    char *readSeq = rsDataSeq(tmpRs);
    // TODO For now, we ignore those unmapped reads.
    if (readRname == NULL) {
      tmpRs = gsItNextRec(gsIt);
      if (tmpRs == NULL) {
        tmpCs = gsItNextChrom(gsIt);
        tmpRs = gsItNextRec(gsIt);
      }
      continue;
    }

    // TODO test the usage of aux_xxx methods from htslib/sam.h
    // TODO see "sam.c" in htslib/test. And check the method "bam_aux_get()"
    // TODO see "SAMtags.pdf" and check the value types when extracting tags.
    // TODO extract tags and print out
    // TODO try modifying the tags and modify the values of those tags

    // --------------------- get the ref sequence -----------------------
    ChromFa *tmpCf = getChromFromGenomeFabyName(readRname, gf);
    // Theses XXXPOS are all 1-based.
    int64_t refStartPos = readStartPos - ecLen;
    int64_t refEndPos = readStartPos + readLength - 1 + ecLen;
    if (refStartPos < 1) refStartPos = 1;
    if (refEndPos > tmpCf->length) refEndPos = tmpCf->length;
    char *refSeq = getSeqFromChromFa(refStartPos, refEndPos, tmpCf);
    // printf("recSam - startPos: %" PRId64 ", endPos: %" PRId64 ", rname:
    // %s\n",
    //        refStartPos, refEndPos, readRname);

    // --------- get all variants' combinations within the interval -----
    // find the first variant that can be integrated and start iteration on
    // ChromVcf to find the rest variants
    tmpCv = getChromFromGenomeVcfbyName(readRname, gv);
    RecVcf *firstRv =
        findFirstVarThatCanBeIntegrated(tmpCv, refStartPos, refEndPos);
    tmpRv = firstRv;

    int integratedRvCnt = 0;
    // 1st loop - calculate number of vcf records that needs integration
    // printf("-----------------1st loop-----------------\n");
    while (tmpRv != NULL) {
      if (countIntegratedAllele(tmpRv, refStartPos, refEndPos) > 0) {
        integratedRvCnt++;
        // printVcfRecord_brief(gv, rvData(tmpRv));
      } else {
        break;
      }
      tmpRv = tmpRv->next;
    }
    // printf("integratedRvCnt: %d\n", integratedRvCnt);

    // 2nd loop - save pointers to vcf records that needs integration, calculate
    // number of alleles of each vcf record that needs integration and get ready
    // for the next step (calculation of combinations)
    // printf("-----------------2nd loop-----------------\n");
    ElementRecVcf **ervArray =
        (ElementRecVcf **)calloc(integratedRvCnt, sizeof(ElementRecVcf *));
    tmpRv = firstRv;
    int integratedRv_idx = 0;
    while (tmpRv != NULL) {
      int integratedAlleleCnt =
          countIntegratedAllele(tmpRv, refStartPos, refEndPos);
      if (integratedAlleleCnt > 0) {
        // printVcfRecord_brief(gv, rvData(tmpRv));
        // If the vcf record can be integrated, find all alleles that needs
        // integration on this record.
        ervArray[integratedRv_idx] =
            (ElementRecVcf *)malloc(sizeof(ElementRecVcf));
        ervArray[integratedRv_idx]->rv = tmpRv;
        ervArray[integratedRv_idx]->alleleIdx =
            (int *)calloc(integratedAlleleCnt, sizeof(int));
        ervArray[integratedRv_idx]->alleleCnt = 0;

        int tmp_alleleIdx = 0;
        for (int i = 0; i < rvDataAlleleCnt(tmpRv); i++) {
          if (ifCanIntegrateAllele(tmpRv, i, refStartPos, refEndPos) == 1) {
            ervArray[integratedRv_idx]->alleleIdx[tmp_alleleIdx] = i;
            ervArray[integratedRv_idx]->alleleCnt++;
            tmp_alleleIdx++;
          }
        }
        integratedRv_idx++;
      } else {
        break;
      }
      tmpRv = tmpRv->next;
    }
    // printf("rvArray: \n");
    // for (int i = 0; i < integratedRvCnt; i++) {
    //   printVcfRecord_brief(gv, rvData(ervArray[i]->rv));
    //   printf("\tintegrated alleles' idxes: ");
    //   for (int j = 0; j < ervArray[i]->alleleCnt; j++) {
    //     printf("%d ", ervArray[i]->alleleIdx[j]);
    //   }
    //   printf("\n");
    // }

    // ---------------------- perform integration -----------------------
    int *ervIdxes = (int *)calloc(integratedRvCnt, sizeof(int));
    for (int i = 0; i < integratedRvCnt; i++) ervIdxes[i] = i;
    int newRecCnt = 0;
    for (int i = 1; i < integratedRvCnt + 1; i++) {
      Combinations *cbs = combinations(ervIdxes, integratedRvCnt, i);
      // printf("erv combi(size: %d) cnt: %d\n", i, cbs->combiCnt);
      for (int j = 0; j < cbs->combiCnt; j++) {
        // printf("erv combi[%d]: ", j);
        // for (int k = 0; k < cbs->combiSize; k++) {
        //   printf("%d ", cbs->combis[j][k]);
        // }
        // printf("\n");
        AlleleCombinations *acbs = alleleCombinations(
            ervArray, integratedRvCnt, cbs->combis[j], cbs->combiSize);
        // printf("allele combi cnt: %d\n", acbs->combiCnt);
        // print_AlleleCombinations(acbs);
        for (int k = 0; k < acbs->combiCnt; k++) {
          newRecCnt++;
          integrateVarAndRealign_refactored(
              ervArray, acbs->rvCombi, acbs->alleleCombis[k], acbs->combiSize,
              refStartPos, refSeq, readSeq, gv, tmpRs, op_file, gs);
        }
        destroy_AlleleCombinations(acbs);
      }
    }
    // printf("new rec cnt: %d\n", newRecCnt);
    // printf("\n");

    // --------------------- free allocated memory ----------------------
    free(refSeq);
    free(readSeq);
    for (int i = 0; i < integratedRvCnt; i++) {
      free(ervArray[i]);
    }
    free(ervArray);
    free(ervIdxes);

    // ----------------------- keep on iterating --------------------
    tmpRs = gsItNextRec(gsIt);
    if (tmpRs == NULL) {
      tmpCs = gsItNextChrom(gsIt);
      tmpRs = gsItNextRec(gsIt);
    }
  }

  // Free structures
  destroy_GenomeSamIterator(gsIt);
  destroy_GenomeVcfIterator(gvIt);

  sam_close(op_file);

  destroy_GenomeFa(gf);
  destroy_GenomeSam(gs);
  destroy_GenomeVcf(gv);

  printf("...successfully finished.\n");
}

// *************************************************************
// *************************************************************
// ******************** other operatins ************************
// *************************************************************
// *************************************************************

void selectBadReads(Options *opts) {
  if (getSamFile(opts) == NULL || getOutputFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments not complete for \'selectBadReads\' option.\n");
    exit(EXIT_FAILURE);
  }
  htsFile *samFile = hts_open(getSamFile(opts), "r");
  htsFile *outputFile = hts_open(getOutputFile(opts), "w");

  sam_hdr_t *samHeader = sam_hdr_read(samFile);
  sam_hdr_t *outputHeader = sam_hdr_init();

  bam1_t *record = bam_init1();

  int initError = 0;
  if (samHeader == NULL || outputHeader == NULL) {
    fprintf(stderr, "Error: failed creating BAM header struct.\n");
    initError = 1;
  }
  if (record == NULL) {
    fprintf(stderr, "Error: out of memory allocating BAM struct.\n");
    initError = 1;
  }
  if (initError) {
    bam_destroy1(record);
    sam_hdr_destroy(samHeader);
    sam_hdr_destroy(outputHeader);
    hts_close(outputFile);
    hts_close(samFile);
    exit(EXIT_FAILURE);
  }

  if (sam_hdr_write(outputFile, samHeader) < 0) exit(EXIT_FAILURE);
  const int threshold = MAPQ_threshold(opts);
  for (int ret = sam_read1(samFile, samHeader, record); ret >= 0;
       ret = sam_read1(samFile, samHeader, record)) {
    // I did not find any method or macros to access "qual", so I directly
    // access using pointers and structures ...
    uint8_t quality = record->core.qual;
    if (quality < threshold)
      if (sam_write1(outputFile, samHeader, record) < 0) exit(EXIT_FAILURE);
  }

  bam_destroy1(record);

  sam_hdr_destroy(samHeader);
  sam_hdr_destroy(outputHeader);

  hts_close(outputFile);
  hts_close(samFile);
}

static inline void check_param_integrateVcfToSam(Options *opts) {
  if (getSamFile(opts) == NULL) {
    fprintf(stderr, "Error: lack input *.sam file. \n");
    exit(EXIT_FAILURE);
  }
  if (getVcfFile(opts) == NULL) {
    fprintf(stderr, "Error: lack iput *.vcf file. \n");
    exit(EXIT_FAILURE);
  }
  if (getFaFile(opts) == NULL) {
    fprintf(stderr, "Error: lack input *.fa file. \n");
    exit(EXIT_FAILURE);
  }
  if (getOutputFile(opts) == NULL) {
    fprintf(stderr, "Error: output file undesignateed. \n");
    exit(EXIT_FAILURE);
  }
}

void integrateVcfToSam(Options *opts) {
  check_param_integrateVcfToSam(opts);

  GenomeFa *gf = init_GenomeFa();
  GenomeSam *gs = init_GenomeSam();
  GenomeVcf *gv = init_GenomeVcf();

  loadGenomeFaFromFile(gf, getFaFile(opts));
  loadGenomeSamFromFile(gs, getSamFile(opts));
  loadGenomeVcfFromFile(gv, getVcfFile(opts));

  // Write header for output file (new sam file)
  samFile *out_file = NULL;
  out_file = sam_open(getOutputFile(opts), "w");
  if (out_file == NULL) {
    fprintf(stderr, "Error: cannot open file %s with mode \"w\"\n",
            getOutputFile(opts));
    exit(EXIT_FAILURE);
  }
  sam_hdr_t *header = gsDataHdr(gs);
  if (sam_hdr_write(out_file, header) < 0) {
    fprintf(stderr, "Error: failed to write sam file header. \n");
    exit(EXIT_FAILURE);
  }

  // Build iterator for sam records and vcf records
  GenomeSamIterator *gsIt = init_GenomeSamIterator(gs);
  ChromSam *tmpCs = gsItNextChrom(gsIt);
  RecSam *tmpRs = gsItNextRec(gsIt);

  GenomeVcfIterator *gvIt = init_GenomeVcfIterator(gv);
  ChromVcf *tmpCv = NULL;
  RecVcf *tmpRv = NULL;

  // Iterate all sam records and process them one by one.
  while (tmpRs != NULL) {
    // ------------- get information of temporary sam record ------------
    const char *readRname = rsDataRname(gs, tmpRs);
    const char *readQname = bam_get_qname(rsData(tmpRs));
    int64_t readStartPos = rsDataPos(tmpRs);
    uint32_t readLength = rsDataSeqLength(tmpRs);
    char *readSeq = rsDataSeq(tmpRs);
    // In this project, we ignore those unmapped reads.
    if (readRname == NULL) {
      tmpRs = gsItNextRec(gsIt);
      if (tmpRs == NULL) {
        tmpCs = gsItNextChrom(gsIt);
        tmpRs = gsItNextRec(gsIt);
      }
      continue;
    }

    // TODO finish your design before implementation

    // ----------------------- keep on iterating --------------------
    tmpRs = gsItNextRec(gsIt);
    if (tmpRs == NULL) {
      tmpCs = gsItNextChrom(gsIt);
      tmpRs = gsItNextRec(gsIt);
    }
  }

  // Free structures
  destroy_GenomeSamIterator(gsIt);
  destroy_GenomeVcfIterator(gvIt);

  sam_close(out_file);

  destroy_GenomeFa(gf);
  destroy_GenomeSam(gs);
  destroy_GenomeVcf(gv);

  printf("...successfully finished.\n");
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/************************* Debug Methods ************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/

static int _test_realignment() {
  const char *tseq = "AAAAAAAAACGTACGTACGTAAAAACCCCCGTGTGA";
  const int tlen = strlen(tseq);
  const int64_t tpos = 101;  // 1-based, included
  const char *qseq = "TTTTACGTACGTACCCCCGTAAA";
  const int qlen = strlen(qseq);
  const int64_t qpos = 112;  // 1-based, included

  const char *old_cigar_str = "3I9M4D8M3I";

  // Transfer the cigar string into coded integer array
  uint32_t *old_cigar_array = NULL;
  char *end = NULL;
  size_t m = 0;
  int old_cigar_len =
      sam_parse_cigar(old_cigar_str, &end, &old_cigar_array, &m);

  // for(int i = 0; i < old_cigar_len; i++){
  //   int opLen = bam_cigar_oplen(old_cigar_array[i]);
  //   char opChr = bam_cigar_opchr(old_cigar_array[i]);
  //   printf("%d%c", opLen, opChr);
  // }
  // printf("\n");

  // TODO aborted test

  return 1;
}

void _testSet_grbvOperations() { return; }
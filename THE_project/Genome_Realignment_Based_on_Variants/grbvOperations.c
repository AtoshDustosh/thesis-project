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

// TODO --------------------------- split line ---------------------------
// TODO --------------------------- split line ---------------------------
// TODO --------------------------- split line ---------------------------
// TODO --------------------------- split line ---------------------------
// TODO --------------------------- split line ---------------------------

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
        // TODO ignore other kinds of varaints
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
        // TODO ignore other kinds of varaints
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
      // TODO ignore other kinds of varaints
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

/**
 * @brief Check whether or not to integrate an SNP.
 *
 * @param varPos 1-based variant position
 * @param refStartPos 1-based position of the start point of refSeq
 * @param refEndPos 1-based position of the end point of refSeq
 */
static inline int ifShouldIntegrateSNP(int64_t varPos, int64_t refStartPos,
                                       int64_t refEndPos) {
  if (varPos >= refStartPos && varPos <= refEndPos) {
    /*
     * ref: ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx---------
     * del: -----------------------------------S----------------------------
     */
    return 1;
  } else {
    return 0;
  }
}

static inline int ifShouldIntegrateINS(int64_t varPos, int64_t refStartPos,
                                       int64_t refEndPos) {
  if (varPos >= refStartPos && varPos <= refEndPos) {
    /*
     * ref: ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx---------
     * del: -----------------------------------IIIIIIIII--------------------
     */
    return 1;
  } else {
    return 0;
  }
}

static inline int ifShouldIntegrateDEL(int64_t varStartPos, int64_t varEndPos,
                                       int64_t refStartPos, int64_t refEndPos) {
  // this is based on the assumption that varEndPos > varStartPos
  if (varEndPos >= refStartPos && varStartPos <= refEndPos) {
    /*
     * ref: ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx---------
     * del: -----------------------------------DDDDDDDDD--------------------
     * del: ---------------------------------------------------DDDDDDDD-----
     * del: ---------DDDDDDD------------------------------------------------
     */
    return 1;
  } else {
    /*
     * ref: ------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx---------
     * del: ---DDDDDDD------------------------------------------------------
     * del: ---------------------------------------------------------DDDD---
     */
    return 0;
  }
}

/**
 * @brief Integrate an allele into the reference genome if there exists common
 * region between the variant and reference genome, and then do the realignment.
 * Do nothing if no common region exists.
 * @return int number of sam records generated after this variant is integrated
 * (>= 0); -1 if error occurs.
 */
static inline int integrateVarAndRealign(RecVcf *rv, RecSam *rs, char *refSeq,
                                         char *readSeq, int64_t refStartPos,
                                         int64_t refEndPos, GenomeVcf *gv) {
  assert(rv != NULL && refSeq != NULL);
  int generatedCnt = 0;
  int alleleCnt = rvData(rv)->n_allele;
  if (alleleCnt == 1) {  // this vcf record contains no ALT but only REF
    return generatedCnt;
  } else {  // this vcf record contains ALT (variants)
    printVcfRecord_brief(gv, rvData(rv));
    int64_t varPos = rvDataPos(rv);
    for (int i = 1; i < alleleCnt; i++) {
      AlignResult *ar = init_AlignResult();
      switch (bcf_get_variant_type(rvData(rv), i)) {
        case VCF_SNP: {
          if (ifShouldIntegrateSNP(varPos, refStartPos, refEndPos) ==
              0) {  // check SNP
            break;
          } else {
            int varIdx = varPos - refStartPos;
            char snp = rvData(rv)->d.allele[i][0];
            char *newSeq = strdup(refSeq);
            newSeq[varIdx] = snp;
            // ---------------------printing------------------------
            printf("old ref: %s\n", refSeq);
            align_ksw2(newSeq, strlen(newSeq), readSeq, strlen(readSeq), ar);
            printf("new ref: %s\n", newSeq);
            printf("od score: %" PRIu8 ", old cigar: ", rsDataMapQ(rs));
            for (int i = 0; i < rsData(rs)->core.n_cigar; i++) {
              uint32_t cigarOpt = bam_get_cigar(rsData(rs))[i];
              uint32_t oplen = bam_cigar_oplen(cigarOpt);
              char opchar = bam_cigar_opchr(cigarOpt);
              printf("%" PRId32 "%c", oplen, opchar);
            }
            printf(", score: %" PRIu8 ", cigar: %s\n", ar->score, ar->cigar);
            // ---------------------free-----------------------------
            free(newSeq);
            generatedCnt++;
            break;
          }
        }
        case VCF_INDEL: {
          if (strlen(rvData(rv)->d.allele[0]) == 1) {  // this is an insertion
            if (ifShouldIntegrateINS(varPos, refStartPos, refEndPos) == 0) {
              break;  // check INS
            } else {
              int varIdx = varPos - refStartPos;
              char *inserted = rvData(rv)->d.allele[i];
              inserted++;  // pass the first base (the REF field)
              char *newSeq = insertStr(refSeq, inserted, varIdx);
              // ---------------------printing------------------------
              printf("old ref: %s\n", refSeq);
              printf("new ref: %s\n", newSeq);
              printf("od score: %" PRIu8 ", old cigar: ", rsDataMapQ(rs));
              for (int i = 0; i < rsData(rs)->core.n_cigar; i++) {
                uint32_t cigarOpt = bam_get_cigar(rsData(rs))[i];
                uint32_t oplen = bam_cigar_oplen(cigarOpt);
                char opchar = bam_cigar_opchr(cigarOpt);
                printf("%" PRId32 "%c", oplen, opchar);
              }
              align_ksw2(newSeq, strlen(newSeq), readSeq, strlen(readSeq), ar);
              printf(", score: %" PRIu8 ", cigar: %s\n", ar->score, ar->cigar);
              // ---------------------free-----------------------------
              free(newSeq);
              generatedCnt++;
            }
          } else {  // this is a deletion
            int64_t varEndPos = varPos + rvDataMaxVarLength(rv) - 1;
            if (ifShouldIntegrateDEL(varPos, varEndPos, refStartPos,
                                     refEndPos) == 0) {
              break;  // check DEL
            } else {
              char *deleted = rvData(rv)->d.allele[0];
              deleted++;
              uint32_t lengthDeleted = strlen(deleted);
              char *newSeq =
                  (char *)calloc(strlen(refSeq) - lengthDeleted, sizeof(char));
              // *Idx is the 0-based index in the refSeq for bases
              int varStartIdx =
                  varPos - refStartPos + 1;  // "+1": ignore the REF
              if (varStartIdx < 0) varStartIdx = 0;
              int varEndIdx = varPos + lengthDeleted - refStartPos;
              if (varEndPos > refEndPos) varEndIdx = strlen(refSeq) - 1;
              for (int i = 0; i < varStartIdx; i++) {
                newSeq[i] = refSeq[i];
              }
              for (int i = varEndIdx; i < strlen(refSeq); i++) {
                newSeq[i - lengthDeleted] = refSeq[i];
              }
              // ---------------------printing------------------------
              printf("old ref: %s\n", refSeq);
              printf("new ref: %s\n", newSeq);
              printf("od score: %" PRIu8 ", old cigar: ", rsDataMapQ(rs));
              for (int i = 0; i < rsData(rs)->core.n_cigar; i++) {
                uint32_t cigarOpt = bam_get_cigar(rsData(rs))[i];
                uint32_t oplen = bam_cigar_oplen(cigarOpt);
                char opchar = bam_cigar_opchr(cigarOpt);
                printf("%" PRId32 "%c", oplen, opchar);
              }
              align_ksw2(newSeq, strlen(newSeq), readSeq, strlen(readSeq), ar);
              printf(", score: %" PRIu8 ", cigar: %s\n", ar->score, ar->cigar);
              // ---------------------free-----------------------------
              free(newSeq);
              generatedCnt++;
            }
          }
          break;
        }
        default: {
          fprintf(stderr,
                  "Warning: integration for variant type not supported. \n");
        }
      }
      destroy_AlignResult(ar);
    }
  }
  return generatedCnt;
}

void integrateVcfToSam(Options *opts) {
  const uint8_t ecLen = 5;  // error control length
  if (opts->samFile == NULL || opts->vcfFile == NULL ||
      opts->outputFile == NULL) {
    fprintf(
        stderr,
        "Error: arguments not complete for \"integrateVcfToSam\" option. \n");
    exit(EXIT_FAILURE);
  }

  GenomeFa *gf = init_GenomeFa();
  GenomeSam *gs = init_GenomeSam();
  GenomeVcf *gv = init_GenomeVcf();

  loadGenomeFaFromFile(gf, opts->faFile);
  loadGenomeSamFromFile(gs, opts->samFile);
  loadGenomeVcfFromFile(gv, opts->vcfFile);

  // printGenomeFa(gf);
  // printGenomeSam(gs);
  // printGenomeVcf(gv);

  // iterate sam and vcf at the same time.
  GenomeSamIterator *gsIt = init_GenomeSamIterator(gs);
  ChromSam *tmpCs = gsItNextChrom(gsIt);
  RecSam *tmpRs = gsItNextRec(gsIt);
  GenomeVcfIterator *gvIt = init_GenomeVcfIterator(gv);
  ChromVcf *tmpCv = NULL;
  RecVcf *tmpRv = NULL;
  uint32_t ifSameCv = 0;

  // Iterate all sam records and process variants
  while (tmpRs != NULL) {
    // -------------- extract data from sam record ------------
    const char *readRname = rsDataRname(gs, tmpRs);
    if (readRname == NULL) {
      tmpRs = gsItNextRec(gsIt);
      if (tmpRs == NULL) {
        tmpCs = gsItNextChrom(gsIt);
        tmpRs = gsItNextRec(gsIt);
      }
      continue;
    }
    char *readQname = bam_get_qname(rsData(tmpRs));
    int64_t readStartPos = rsDataPos(tmpRs);
    uint32_t readLength = rsDataSeqLength(tmpRs);
    char *readSeq = rsDataSeq(tmpRs);

    // --------------- get the ref sequence ------------------
    ChromFa *tmpCf = getChromFromGenomeFabyName(readRname, gf);
    int64_t refStartPos = readStartPos - ecLen;
    int64_t refEndPos = readStartPos + readLength - 1 + ecLen;
    if (refStartPos < 1) refStartPos = 1;
    if (refEndPos > tmpCf->length) refEndPos = tmpCf->length;
    char *refSeq = getSeqFromChromFa(refStartPos, refEndPos, tmpCf);

    printf("readQName: %s, readStartPos: %" PRIu64 ", readLen: %" PRIu32
           ", readRname: %s\n",
           readQname, readStartPos, readLength, readRname);
    printf("refStartPos: %" PRId64 ", refEndPos: %" PRId64 "\n", refStartPos,
           refEndPos);

    // -------------- locate the valid variants --------------
    // get the ChromVcf
    ifSameCv = 0;
    if (tmpCv != NULL) {
      if (strcmp(readRname, tmpCv->name) != 0) {
        tmpCv = getChromFromGenomeVcfbyName(readRname, gv);
        if (tmpCv == NULL) {
          fprintf(stderr, "Error: didn't found chrom with name: %s\n",
                  readRname);
          exit(EXIT_FAILURE);
        }
      } else {
        ifSameCv = 1;
      }
    } else {
      tmpCv = getChromFromGenomeVcfbyName(readRname, gv);
    }

    // 1. find the first variant that should be integrated
    // 2. integrate the variant and check whether the next vaiant should be
    // integrated
    tmpRv = findFirstVarThatCanBeIntegrated(tmpCv, refStartPos, refEndPos);
    while (tmpRv != NULL &&
           ifCanIntegrateVar(tmpRv, refStartPos, refEndPos) == 1) {
      if (integrateVarAndRealign(tmpRv, tmpRs, refSeq, readSeq, refStartPos,
                                 refEndPos, gv) < 0) {
        fprintf(stderr, "Warning: integration failed for this record.\n");
        printVcfRecord_brief(gv, rvData(tmpRv));
        printf("refStartPos: %" PRId64 ", refEndPos: %" PRId64 "\n",
               refStartPos, refEndPos);
      }
      tmpRv = tmpRv->next;
    }
    printf("\n");

    // -------------------------- split line -------------------
    free(readSeq);
    free(refSeq);

    // the following codes are the core of sam-records-iteration
    tmpRs = gsItNextRec(gsIt);
    if (tmpRs == NULL) {
      tmpCs = gsItNextChrom(gsIt);
      tmpRs = gsItNextRec(gsIt);
    }
  }

  destroy_GenomeSamIterator(gsIt);
  destroy_GenomeVcfIterator(gvIt);

  destroy_GenomeFa(gf);
  destroy_GenomeSam(gs);
  destroy_GenomeVcf(gv);
}

// TODO refactoring ...
void integrateVcfToSam_refactored(Options *opts) {
  if (opts->samFile == NULL || opts->vcfFile == NULL) {
    fprintf(stderr,
            "Error: arguments not complete for variants integration. \n");
    exit(EXIT_FAILURE);
  }
  if (opts->outputFile == NULL) {
    opts->outputFile = "data/defaultOutput.txt";
  }

  GenomeFa *gf = init_GenomeFa();
  GenomeSam *gs = init_GenomeSam();
  GenomeVcf *gv = init_GenomeVcf();

  loadGenomeFaFromFile(gf, opts->faFile);
  loadGenomeSamFromFile(gs, opts->samFile);
  loadGenomeVcfFromFile(gv, opts->vcfFile);

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

    // --------------------- get the ref sequence -----------------------
    ChromFa *tmpCf = getChromFromGenomeFabyName(readRname, gf);
    int64_t refStartPos = readStartPos - ecLen;
    int64_t refEndPos = readStartPos + readLength - 1 + ecLen;
    if (refStartPos < 1) refStartPos = 1;
    if (refEndPos > tmpCf->length) refEndPos = tmpCf->length;
    char *refSeq = getSeqFromChromFa(refStartPos, refEndPos, tmpCf);

    // TODO
    // --------- get all variants' combinations within the interval -----
    // find the first variant that can be integrated and start iteration on
    // ChromVcf to find the rest variants
    tmpCv = getChromFromGenomeVcfbyName(readRname, gv);
    RecVcf *firstRv =
        findFirstVarThatCanBeIntegrated(tmpCv, refStartPos, refEndPos);
    tmpRv = firstRv;

    int integratedRvCnt = 0;
    int *integratedAlleleCnts = NULL;
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
    // TODO temporary  task:
    // TODO print the results for checking correctness
    // TODO temporary  task:
    // TODO print the results for checking correctness
    // TODO temporary  task:
    // TODO print the results for checking correctness
    // TODO temporary  task:
    // TODO print the results for checking correctness
    // TODO temporary  task:
    // TODO print the results for checking correctness


    // 3rd loop - create PendingAlleles for next step of calculating
    // combinations
    // printf("-----------------3rd loop-----------------\n");

    printf("integratedRvCnt: %d\n", integratedRvCnt);
    printf("recSam - startPos: %" PRId64 ", endPos: %" PRId64 ", rname: %s\n",
           refStartPos, refEndPos, readRname);
    printf("rvArray: \n");
    for (int i = 0; i < integratedRvCnt; i++) {
      // printVcfRecord_brief(gv, rvData(rvArray[i]));
    }
    // for (int i = 0; i < integratedRvCnt; i++) {
    //   PendingAlleles *tmpAllele = alleles_pending[i];
    //   printf("pending allele (%d) idxes: (%d) ", i, tmpAllele->rvIdx);
    //   for (int j = 0; j < tmpAllele->alleleCnt; j++) {
    //     printf("%d ", tmpAllele->alleleIdx[j]);
    //     // TODO this part of code is not correct
    //     // printing of "PendingAlleles" unfinished
    //   }
    //   printf("\n");
    // }

    // -------------- generate all permutaions of variants --------------

    // --------------------- perform integration ------------------------

    // ----------------------- keep on iterating --------------------
    printf("\n");
    tmpRs = gsItNextRec(gsIt);
    if (tmpRs == NULL) {
      tmpCs = gsItNextChrom(gsIt);
      tmpRs = gsItNextRec(gsIt);
    }
  }

  destroy_GenomeFa(gf);
  destroy_GenomeSam(gs);
  destroy_GenomeVcf(gv);
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
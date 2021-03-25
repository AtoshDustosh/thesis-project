#include "grbvOperations.h"

// TODO Static methods for integration. Consider encapsulating ...

/**
 * @brief Judge whether a variant should be integrated or not.
 *
 * @param rv variant
 * @param startPos 1-based position of the start point of refSeq
 * @param endPos 1-based position of the end point of refSeq
 * @return int 1 if the variant should be integrated; 0 otherwise.
 */
static inline int ifShouldIntegrateVar(RecVcf *rv, int64_t startPos,
                                       int64_t endPos) {
  // TODO
  int64_t varStartPos = rvDataPos(rv);
  int64_t varEndPos = varStartPos + rvDataMaxVarLength(rv);
  if (varStartPos < endPos && varEndPos > startPos)
    return 1;
  else
    return 0;
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
static inline RecVcf *findFirstVarToIntegrate(ChromVcf *cv, int64_t startPos,
                                              int64_t endPos) {
  assert(cv != NULL);
  ChromVcfIterator *cvIt = init_ChromVcfIterator(cv);
  RecVcf *rv = cvItNextRec(cvIt);
  while (rv != NULL) {
    if (ifShouldIntegrateVar(rv, startPos, endPos)) {
      destroy_ChromVcfIterator(cvIt);
      return rv;
    }
    rv = cvItNextRec(cvIt);
  }

  destroy_ChromVcfIterator(cvIt);
  return NULL;
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
  /*
   * ref: --------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx-------
   * read: -------------|-|-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx------|-------
   * refStartPos: ------#-|---------------------------------------------|-------
   * refEndPos: ----------|---------------------------------------------#-------
   * example SNP: --------Y-----------------------------------------------------
   * varPos - refStartPos = 2
   * idx in the refSeq = 2
   * modify refSeq[2] from "X" to "Y"
   */
  assert(rv != NULL && refSeq != NULL);
  int generatedCnt = 0;
  int alleleCnt = rvData(rv)->n_allele;
  if (alleleCnt == 1) { // this vcf record contains no ALT but only REF
    return generatedCnt;
  } else {  // this vcf record contains ALT (variants)
    // printf("integrated: ");
    // printVcfRecord_brief(gv, rvData(rv));
    static char newSeq[MAX_INFO_LENGTH];
    int64_t varPos = rvDataPos(rv);
    for (int i = 1; i < alleleCnt; i++) {
      switch (bcf_get_variant_type(rvData(rv), i)) {
        case VCF_SNP: {
          // TODO
          int varIdx = varPos - refStartPos;
          char snp = rvData(rv)->d.allele[i][0];
          break;
        }
        case VCF_INDEL: {
          // TODO
          break;
        }
        default: {
          fprintf(stderr,
                  "Warning: integration for variant type not supported. \n");
        }
      }
    }
  }
  return generatedCnt;
}

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

  // TODO iterate sam and vcf at the same time.
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
    printf("readSeq: %s\n", readSeq);
    printf("refSeq:  %s\n", refSeq);
    printf("\n");

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
    tmpRv = findFirstVarToIntegrate(tmpCv, refStartPos, refEndPos);
    while (tmpRv != NULL &&
           ifShouldIntegrateVar(tmpRv, refStartPos, refEndPos)) {
      if (integrateVarAndRealign(tmpRv, tmpRs, refSeq, readSeq, refStartPos,
                                 refEndPos, gv) < 0) {
        fprintf(stderr, "Warning: integration failed for this record.\n");
        printVcfRecord_brief(gv, rvData(tmpRv));
        printf("refStartPos: %" PRId64 ", refEndPos: %" PRId64 "\n",
               refStartPos, refEndPos);
      }
      tmpRv = tmpRv->next;
    }

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
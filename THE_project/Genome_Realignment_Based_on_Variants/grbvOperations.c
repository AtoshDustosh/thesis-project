#include "grbvOperations.h"


// TODO Static methods for integration. Consider encapsulating ...
/**
 * @brief Judge whether should integrate a variant onto the reference genome. 
 * 
 * @param rv variant
 * @param startPos 1-based position of integrated area's start point
 * @param endPos 1-based position of integrated area's end point
 * @return int If should, return 1; 0 otherwise. 
 */
static inline int shouldIntegrateVar(RecVcf *rv, int64_t startPos, int64_t endPos){
  // TODO
  return 0;
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
  const uint8_t ecLen = 15;  // error control length
  if (opts->samFile == NULL || opts->vcfFile == NULL) {
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

  // Iterate all sam records
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
        tmpCv = getChromFromGenomeVcf(readRname, gv);
      } else {
        ifSameCv = 1;
      }
    } else {
      tmpCv = getChromFromGenomeVcf(readRname, gv);
    }
    
    // TODO the implementation below is too redundant and buggy. Considering new method for implementation ...
    // find variants that have common bases with the reads.
    // TODO the following codes are abandoned. Trying new implementation. 
    RecVcf *startRv = getRecBeforePosFromChromVcf(readStartPos, tmpCv);
    RecVcf *endRv =
        getRecAfterPosFromChromVcf(readStartPos + readLength, tmpCv);
    while (startRv != endRv) {
      // TODO Multiple situations. Could be the following cases:
      // Ref: --------XXXXXXXXXXXXXXXXXXXXXXXXXXXX----------
      // Read:-----------YYYYYYYYYYYYYYYYYYYYY--------------
      // var: -----AAAAAAAAA--------------------------------
      // var: ----------------BBBBBBBB----------------------
      // var: ----------------------------CCCCCCCCCC--------
      // var: ----DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD---
      uint64_t rvStartPos = rvDataPos(startRv);
      uint64_t rvEndPos = rvStartPos + rvDataMaxVarLength(startRv) - 1;
      if (rvStartPos > refStartPos && rvEndPos < refEndPos) {
        // TODO
      } else if (rvStartPos < refStartPos && rvEndPos > refStartPos) {
        // TODO
      } else if (rvStartPos > refStartPos && rvEndPos > refEndPos) {
        // TODO
      } else if (rvStartPos < refStartPos && rvEndPos > refEndPos) {
        // TODO
      } else {
        // actually there are many cases not listed here. But we just consider
        // the ordinary cases for now ... temporarily
      }
    }

    // -------------------------- split line -------------------
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
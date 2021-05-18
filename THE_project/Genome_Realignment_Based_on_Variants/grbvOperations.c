#include "grbvOperations.h"

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
#include "grbvOperations.h"

/**
 * @brief Check whether opts is valid. 
 * 
 * @return int 1 if valid; 0 otherwise
 */
static int checkOpt_selectBadReads(Options *opts)
{
  if (getSamFile(opts) == NULL || getOutputFile(opts) == NULL)
  {
    fprintf(stderr, "Error: arguments not complete for \'selectBadReads\' option.\n");
    return 0;
  }
  return 1;
}

void selectBadReads(Options *opts)
{
  if (!checkOpt_selectBadReads(opts))
    return;
  htsFile *samFile = hts_open(getSamFile(opts), "r");
  htsFile *outputFile = hts_open(getOutputFile(opts), "w");

  sam_hdr_t *samHeader = sam_hdr_read(samFile);
  sam_hdr_t *outputHeader = sam_hdr_init();

  bam1_t *record = bam_init1();

  int initError = 0;
  if (samHeader == NULL || outputHeader == NULL)
  {
    fprintf(stderr, "Error: failed creating BAM header struct.\n");
    initError = 1;
  }
  if (record == NULL)
  {
    fprintf(stderr, "Error: out of memory allocating BAM struct.\n");
    initError = 1;
  }
  if (initError)
  {
    bam_destroy1(record);
    sam_hdr_destroy(samHeader);
    sam_hdr_destroy(outputHeader);
    hts_close(outputFile);
    hts_close(samFile);
    exit(EXIT_FAILURE);
  }

  if (sam_hdr_write(outputFile, samHeader) < 0)
    exit(EXIT_FAILURE);
  const int threshold = MAPQ_threshold(opts);
  for (int ret = sam_read1(samFile, samHeader, record); ret >= 0; ret = sam_read1(samFile, samHeader, record))
  {
    // I did not find any method or macros to access "qual", so I directly access using pointers and structures ...
    uint8_t quality = record->core.qual;
    if (quality < threshold)
      if (sam_write1(outputFile, samHeader, record) < 0)
        exit(EXIT_FAILURE);
  }

  bam_destroy1(record);

  sam_hdr_destroy(samHeader);
  sam_hdr_destroy(outputHeader);

  hts_close(outputFile);
  hts_close(samFile);
}
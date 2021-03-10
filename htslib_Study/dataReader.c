#include "dataReader.h"


int vcf_scan_and_print(htsFile *inputFile, Options *opts) {
  bcf_hdr_t *header = bcf_hdr_read(inputFile);
  bcf1_t *record = bcf_init1();

  int initError = 0;
  if (header == NULL) {
    fprintf(stderr, "Failed creating BCF header struct\n");
    initError = 1;
  }
  if (record == NULL) {
    fprintf(stderr, "Out of memory allocating VCF struct\n");
    initError = 1;
  }
  if (initError) {
    bcf_hdr_destroy(header);
    bcf_destroy1(record);
    return EXIT_FAILURE;
  }

  printf("successfully got into the vcf loop\n");

  clock_t start = clock(), end = 0;
  int ret = 0;
  printVcfHeader(header);
  while (ret = bcf_read1(inputFile, header, record) >= 0) {
    /*
     * This function (bcf_unpack()) must be called when reading vcf/bcf file.
     * No need to worry about duplicate work.
     */
    // typedef struct _define_StructTest {
    //   char *name;
    // } StructTest;
    bcf_unpack(record, BCF_UN_ALL);
    printVcfRecord_brief(header, record);
    // StructTest *st = (StructTest *)malloc(sizeof(StructTest));
    // st->name = strdup(bcf_seqname_safe(header, record));
    // printf("chrom name: %s\n", st->name);
    // free(st->name);
    // free(st);
  }
  printf("\n");
  end = clock();
  printf("total time(s): %f\n", (double)(end - start) / CLOCKS_PER_SEC);

  bcf_hdr_destroy(header);
  bcf_destroy1(record);
  return EXIT_SUCCESS;
}

int sam_scan_and_print(htsFile *inputFile, Options *opts) {
  sam_hdr_t *header = sam_hdr_read(inputFile);
  bam1_t *record = bam_init1();

  int initError = 0;
  if (header == NULL) {
    fprintf(stderr, "Failed creating BAM header struct\n");
    initError = 1;
  }
  if (record == NULL) {
    fprintf(stderr, "Out of memory allocating BAM struct\n");
    initError = 1;
  }
  if (initError) {
    sam_hdr_destroy(header);
    bam_destroy1(record);
    return EXIT_FAILURE;
  }

  printf("successfully got into the sam loop\n");

  clock_t start = clock(), end = 0;
  int ret = 0;
  printSamHeader(header);
  while (ret = sam_read1(inputFile, header, record) >= 0) {
    printSamRecord_brief(header, record);
  }
  printf("\n");
  end = clock();
  printf("total time(s): %f\n", (double)(end - start) / CLOCKS_PER_SEC);

  bam_hdr_destroy(header);
  bam_destroy1(record);
  return EXIT_SUCCESS;
}

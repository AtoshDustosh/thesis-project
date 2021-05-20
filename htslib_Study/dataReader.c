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

    // if(record->n_allele == 1){
    //   continue;
    // }else{
    //   if(record->d.allele[1][0] == '<'){
    //     // ignore
    //     printf("ignored %s\n", record->d.allele[1]);
    //     continue;
    //   }
    // }

    // char *dst_str1 = NULL;
    // int ndst = 0;
    // bcf_get_info_string(header, record, "VARIANT_OVERALL_TYPE", &dst_str1,
    //                     &ndst);
    // fprintf(stderr, "VARIANT_OVERALL_TYPE(%d): %s |\t", ret, dst_str1);
    // free(dst_str1);

    // char *dst_str2 = (char *)calloc(10, sizeof(char));
    // ret = bcf_get_info_string(header, record, "SVTYPE", &dst_str2, &ndst);
    // fprintf(stderr, "SVTYPE(%d): %s |\t", ndst, dst_str2);
    // free(dst_str2);

    // int32_t *dst_int32 = (int32_t *)calloc(10, sizeof(int32_t));
    // ret = bcf_get_info_int32(header, record, "SVLEN", &dst_int32, &ndst);
    // fprintf(stderr, "SVLEN(%d): ", ndst);
    // for (int i = 0; i < ndst; i++) {
    //   if (dst_int32[i] == 0) {
    //     break;
    //   } else {
    //     fprintf(stderr, "%" PRId32 ",", dst_int32[i]);
    //   }
    // }
    // printf("|");
    // printf("\n");
    // fprintf(stderr, "%s\t%" PRId64 "\t%s\t",
    //         bcf_seqname_safe(header, record), record->pos + 1, record->d.id);
    // if (record->n_allele == 1) {
    //   fprintf(stderr, ".");
    // } else {
    //   for (int i = 1; i < record->n_allele; i++) {
    //     fprintf(stderr, "%s(%d),", record->d.allele[i],
    //             bcf_get_variant_type(record, i));
    //   }
    // }
    // fprintf(stderr, "(ignore)");
    // fprintf(stderr, "\t%f\t\n", record->qual);
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

    // sam_format_aux1
    char *aux_data = "abcdefg";
    char aux_type = 'Z';
    char *aux_tag = "XV";
    bam_aux_append(record, aux_tag, aux_type, strlen(aux_data) + 1, aux_data);

    printf("modified: ");
    printSamRecord_brief(header, record);
    printf(".......................\n");
  }
  printf("\n");
  end = clock();
  printf("total time(s): %f\n", (double)(end - start) / CLOCKS_PER_SEC);

  bam_hdr_destroy(header);
  bam_destroy1(record);
  return EXIT_SUCCESS;
}

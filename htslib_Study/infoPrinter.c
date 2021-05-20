#include "infoPrinter.h"

void printVcfHeader(bcf_hdr_t *header) {
  printf("vcf version: %s\n", bcf_hdr_get_version(header));
  printf("\n");
}

void printVcfRecord(bcf1_t *record) {
  /*
   * The ids of chroms correspond to the sequence they appear in the
   * vcf/bcf file. For example, if contigs are written as follows,
   * ##contig=<ID=5>
   * ##contig=<ID=1>
   * ##contig=<ID=3>
   * id array of the 3 chroms/contigs should be 0, 1, 2, instead of
   * 3, 1, 2.
   */
  printf("chrom: %d\t", record->rid);
  /*
   * The output on console looks like 0-based index. But that's
   * for the convenience of array in C.
   */
  printf("pos: %ld\t", record->pos);
  printf("qual: %f\n", record->qual);
  /*
   * Including the reference allele, which is "record->d.allele[0]".
   */
  printf("number of alleles: %d\n", record->n_allele);
  printf("alleles: ");
  for (int i = 0; i < record->n_allele; i++) {
    printf("%s ", record->d.allele[i]);
  }
  printf("\n");

  printf("variant types:");
  for (int i = 1; i < record->n_allele; i++) {
    printf("%d ", bcf_get_variant_type(record, i));
  }
  printf("\n");
  // other fields are temporarily not needed
  printf("\n");
}

void printVcfRecord_brief(bcf_hdr_t *hdr, bcf1_t *record) {
  // chrom
  printf("%s\t", bcf_seqname_safe(hdr, record));
  // pos
  printf("%" PRId64 "\t", record->pos);
  // id
  printf("%s\t", record->d.id);
  // ref alt
  printf("%s\t", record->d.allele[0]);
  if (record->n_allele == 1) {
    printf(".");
  } else {
    for (int i = 1; i < record->n_allele; i++) {
      if (i >= 2) {
        printf(",%s", record->d.allele[i]);
      } else {
        printf("%s", record->d.allele[i]);
      }
    }
  }
  printf("\t");
  // qual
  printf("%f\t", record->qual);
  // filter, info, format and other fields are ignored
  printf("\n");
}

void printSamHeader(bam_hdr_t *header) {
  printf("reference count: %d\n", header->ref_count);
  printf("number of reference sequences: %d\n", header->n_targets);
  printf("name and length of reference sequences:\n");
  for (int i = 0; i < header->n_targets; i++) {
    printf("\t%s\t%d\n", header->target_name[i], header->target_len[i]);
  }
  printf("\n");
}

void printSamRecord(bam1_t *record) {
  printf("seq length: %d\n", record->core.l_qseq);
  /*
   * The output on console looks like 0-based index. But that's
   * for the convenience of array in C.
   */
  printf("pos: %ld\n", record->core.pos);
  /*
   * The output is not the same as in the original *.sam file.
   * Chrom IDs are re-ordered by the htslib.
   */
  printf("chrom id: %d\n", record->core.tid);
  printf("flag: 0x%08x\n", record->core.flag);
  printf("sequence length: %d\n", record->core.l_qseq);
  /*
   * It is recommended to use bam_seqi() macro instead of array
   * manipulations on array retrieved by bam_get_seq().
   */
  printf("seq(0x, 1-A, 2-C, 4-G, 8-T, 15-N):\n\t");
  for (int i = 0; i < record->core.l_qseq; i++) {
    printf("%x", bam_seqi(bam_get_seq(record), i));
  }
  printf("\n");
  /*
   * For example, cigar "22S86M" - cigar operation is 2, and
   * the cigar array is a 2-element array.
   * The element is coded in the following format, and each
   * element corresponds to a cigar operation:
   *   (0x) [count of cigar operation][integer operation type code]
   * Thus, "22S86M" may look like "0x 0164 0560"
   */
  printf("cigar operations type count: %d\n", record->core.n_cigar);
  printf("cigar(0x):\n\t");
  for (int i = 0; i < record->core.n_cigar; i++) {
    printf("0x%04x ", bam_get_cigar(record)[i]);
  }
  printf("\n");

  printf("qname: %s\n", bam_get_qname(record));

  /*
   * "+33" is necessary. The original value in the array is as
   * the BAM specification instead of the SAM ASCII printable method.
   */
  printf("qual: \n\t");
  for (int i = 0; i < record->core.l_qseq; i++) {
    printf("%c", bam_get_qual(record)[i] + 33);
  }
  printf("\n");

  printf("\n");
}

void printSamRecord_brief(bam_hdr_t *hdr, bam1_t *record) {
  // qname
  printf("%s\t", bam_get_qname(record));
  // flag
  printf("0x%" PRIx16 "\t", record->core.flag);
  // rname
  printf("%s\t", sam_hdr_tid2name(hdr, record->core.tid));
  // pos
  printf("%" PRId64 "\t", record->core.pos);
  // mapq
  printf("%" PRIu8 "\t", record->core.qual);
  // cigar
  for (int i = 0; i < record->core.n_cigar; i++) {
    uint32_t cigarOpt = bam_get_cigar(record)[i];
    uint32_t oplen = bam_cigar_oplen(cigarOpt);
    char opchar = bam_cigar_opchr(cigarOpt);
    printf("%" PRId32 "%c", oplen, opchar);
  }
  printf("\t");
  // rnext pnext
  printf("%" PRId32 "\t", record->core.mtid);
  // pnext
  printf("%" PRId64 "\t", record->core.mpos);
  // tlen ... this field is not necessary for this project and there has been
  // some differences on the definition of this field
  printf("*\t");
  // seq (long string)
  for (int i = 0; i < record->core.l_qseq; i++) {
    switch (bam_seqi(bam_get_seq(record), i)) {
      case 1: {
        printf("A");
        break;
      }
      case 2: {
        printf("C");
        break;
      }
      case 4: {
        printf("G");
        break;
      }
      case 8: {
        printf("T");
        break;
      }
      case 15: {
        printf("N");
        break;
      }
    }
  }
  printf("\t");
  // qual (long string)
  for (int i = 0; i < record->core.l_qseq; i++) {
    printf("%c", bam_get_qual(record)[i] + 33);
  }
  printf("\t");
  // "XV" aux field
  printf("XV:%s", bam_aux_get(record, "XV"));
  printf("\t");
  printf("\n");
}

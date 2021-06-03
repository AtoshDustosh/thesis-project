#include "kmerGeneration.h"

static int kmerLength = 0;

/**
 * @brief  Check if the kmer matches specification for output. Valid kmers
 * doesn't contain ignored bases like 'N'.
 * @retval true if the kmer doesn't contains ignored bases; false otherwise.
 */
static inline bool check_kmer_valid(char *kmer) {
  while (*kmer != '\0') {
    for (int i = 0; i < cnt_ignored_bases; i++) {
      if (*kmer == ignored_bases[i]) {
        return false;
      }
    }
    kmer++;
  }
  return true;
}

static inline bool ifCanIntegrateAllele(RecVcf_bplus *rv, int alleleIdx,
                                        int32_t startPos, int32_t endPos) {
  if (alleleIdx == 0) return false;  // ignore REF allele
  const char *var_alt = rv_allele(rv, alleleIdx);
  const char *var_ref = rv_allele(rv, 0);
  int length_alt = strlen(var_alt);
  int length_ref = strlen(var_ref);
  int32_t varStartPos = rv_pos(rv);
  int32_t varEndPos = varStartPos + length_ref - 1;
  if (varStartPos <= endPos && varEndPos >= startPos) {
    return true;
  } else {
    return false;
  }
}

static inline int count_integrated_allele(RecVcf_bplus *rv, int32_t startPos,
                                          int32_t endPos) {
  int varStartPos = rv_pos(rv);
  int cnt_integrated_allele = 0;
  for (int i = 0; i < rv_alleleCnt(rv); i++) {
    switch (bcf_get_variant_type(rv_object(rv), i)) {
      case VCF_REF: {
        break;
      }
      default: {
        const char *var_ref = rv_allele(rv, 0);
        const char *var_alt = rv_allele(rv, i);
        int length_ref = strlen(var_ref);
        int length_alt = strlen(var_alt);
        int32_t varEndPos = varStartPos + length_ref - 1;
        if (varEndPos >= startPos && varStartPos <= endPos)
          cnt_integrated_allele++;
        break;
      }
    }  // end of swich
  }    // end of for
  return cnt_integrated_allele;
}

/**
 * @brief  Generate kmers integrated with varaints and store in the hash table.
 * @param  *ervArray[]: array of varaints on this interval
 * @param  ervCombi[]: array of combination of selected variants
 * @param  alleleCombi[]: array of combinations of alleles on the selected vars
 * @param  length_combi: length of combinations
 * @param  lbound: 1-based expanded lbound of the ref sequence (included)
 * @param  rbound: 1-based expanded rbound of the ref sequence (included)
 * @param  pos_start: 1-based start position of the interval
 * @param  pos_end: 1-based end position of the interval
 * @param  *hashTable: kmers' hash table
 * @param  *cf: structure for a chromosome from reference genome
 * @param  *gv: structure for variants
 * @param  *fp_op: FILE pointer to output file
 */
static inline void generateKmers_integrated(
    Element_RecVcf *ervArray[], int ervCombi[], int alleleCombi[],
    int length_combi, int32_t lbound, int32_t rbound, int32_t pos_start,
    int32_t pos_end, KmerHashTable *hashTable, ChromFa *cf, GenomeVcf_bplus *gv,
    FILE *fp_op) {
  // -------------- set bounds of ref sequence and update bounds ---------------
  const int32_t length_chrom = chromFa_length(cf);
  // expand by 1 base in order to get input char and output char
  int lbound_ref = pos_start - (kmerLength - 1) - 1;  // 1-based
  int rbound_ref = pos_end + (kmerLength - 1) + 1;    // 1-based
  lbound_ref = lbound_ref > 0 ? lbound_ref : 1;

  int length_seq_ref = 0;
  for (int i = 0; i < length_combi; i++) {
    RecVcf_bplus *tmp_rv = ervArray[ervCombi[i]]->rv;
    const char *allele_ref = rv_allele(tmp_rv, 0);
    const char *allele_alt = rv_allele(tmp_rv, alleleCombi[i]);
    int length_allele_ref = strlen(allele_ref);
    int length_allele_alt = strlen(allele_alt);
    printf("allele_ref[%d]: %s, allele_alt[%d]: %s\n", i, allele_ref, i,
           allele_alt);
    // genomeVcf_bplus_printRec(gv,tmp_rv);

    if (length_allele_ref == length_allele_alt) {
      if (length_allele_ref == 1) {  // SNP
        continue;
      } else {  // MNP, but length_alt == length_ref
        length_seq_ref += length_allele_ref;
      }
    } else if (length_allele_alt > length_allele_ref) {  // INS
      length_seq_ref += (length_allele_alt - length_allele_ref);
    } else {  // others
      rbound_ref += (length_allele_ref - length_allele_alt);
    }
  }
  rbound_ref = rbound_ref <= length_chrom ? rbound_ref : length_chrom;
  length_seq_ref += rbound_ref - lbound_ref + 1;
  if (length_seq_ref <= 0) {
    // Do nothing - actually this case should not happen
    fprintf(stderr,
            "Warning: unexpected case happened - length_seq_ref <= 0.\n");
    return;
  }

  // Get original ref sequence
  const char *original_seq_ref = getSeqFromChromFa(lbound_ref, rbound_ref, cf);

  // ------------- extract kmer, input char, output char, and pos --------------
  char buf_kmer[kmerLength + 1];  // buffer for constructing a kmer
  memset(buf_kmer, 0, kmerLength + 1);
  int idx_allele = 0; // index for the allele being handled
  int pos_kmer_ref = 0; // 1-based, pos on the reference sequence
  int pos_kmer_offset = 0;  // for handling kmers within an INS only
  

}

/**
 * @brief  Process and generate kmers for interval [pos_start, pos_end] on the
 * id_chrom-th chromosome. And then output the generated kmers into the
 * specified output file.
 * @param  id_chrom: index of the chromosome
 * @param  pos_start: 1-based start position of interval. included
 * @param  pos_end: 1-based end position of interval. included
 */
static inline void generateKmers_process(int32_t id_chrom, int32_t pos_start,
                                         int32_t pos_end, GenomeFa *gf,
                                         GenomeVcf_bplus *gv, FILE *fp_op) {
  // -------- Expand the interval towards both sides by (kmerLength - 1) -------
  // Fixed interval will ignore kmers that don't have input char or output char
  ChromFa *chrom = getChromFromGenomeFabyIndex(id_chrom, gf);
  int32_t length_chrom = chromFa_length(chrom);
  const char *name_chrom = chromFa_name(chrom);
  assert(length_chrom >= 0);

  // Fix pos to handle kmers at the beginning or end of the chrom
  pos_start = pos_start > 1 ? pos_start : 2;
  pos_end = pos_end < length_chrom ? pos_end : length_chrom - 1;

  // Expand interval (also ignore the kmers at the ends of the chrom)
  int32_t lbound = pos_start - (kmerLength - 1);
  lbound = lbound > 1 ? lbound : 2;
  int32_t rbound = pos_end + (kmerLength - 1);
  rbound = rbound < length_chrom ? rbound : length_chrom - 1;

  // ------------------------ Create hash table for kmers ----------------------
  // Maybe a larger size of the table would be better, but this one can do
  int32_t tableSize = (rbound - lbound);
  KmerHashTable *hashTable = init_kmerHashTable(tableSize, kmerLength);

  // -------------------- Save original kmers into hash table ------------------
  char *seq_ref = getSeqFromChromFa(lbound, rbound, chrom);
  int32_t tmp_offset = 0;
  while (lbound + tmp_offset <= pos_end) {
    int32_t local_lbound = lbound + tmp_offset;
    char inputChar = charOfBase(getBase(chrom, local_lbound - 1));
    char outputChar = charOfBase(getBase(chrom, local_lbound + kmerLength));
    char *kmer = subStr(seq_ref, tmp_offset, tmp_offset + kmerLength - 1);
    if (check_kmer_valid(kmer)) {
      kmerHashTable_add(kmer, lbound + tmp_offset, inputChar, outputChar,
                        hashTable);
      // printf("kmer: %s, pos: %d\n", kmer, (int)(pos_start + tmp_offset));
    }
    tmp_offset++;
  }
  free(seq_ref);

  // ---------- Get variants within the interval [pos_start, pos_end] ----------
  RecVcf_bplus *rv_tmp =
      genomeVcf_bplus_getRecAfterPos(gv, name_chrom, pos_start);
  RecVcf_bplus *first_rv = rv_tmp;

  // 1st loop - calculate number of vcf records that needs integration
  int cnt_integrated_variants = 0;
  while (rv_tmp != NULL) {
    if (count_integrated_allele(rv_tmp, pos_start, pos_end) > 0) {
      cnt_integrated_variants++;
    } else {
      if (rv_pos(rv_tmp) > pos_end) {
        break;
      }
    }
    rv_tmp = next_RecVcf_bplus(rv_tmp);
  }
  // 2nd loop - save pointers to vcf records that needs integration and
  // calculate number of alleles of each vcf record that needs integration.
  int cnt_integrated_alleles = 0;
  Element_RecVcf **ervArray = (Element_RecVcf **)calloc(
      cnt_integrated_variants, sizeof(Element_RecVcf *));
  rv_tmp = first_rv;
  int idx_integrated_variants = 0;
  while (rv_tmp != NULL) {
    int tmp_cnt_integrated_allele =
        count_integrated_allele(rv_tmp, pos_start, pos_end);
    if (tmp_cnt_integrated_allele > 0) {
      // Find all alleles that needs integration on this vcf record
      ervArray[idx_integrated_variants] =
          (Element_RecVcf *)malloc(sizeof(Element_RecVcf));
      ervArray[idx_integrated_variants]->rv = rv_tmp;
      ervArray[idx_integrated_variants]->alleleIdx =
          (int *)calloc(tmp_cnt_integrated_allele, sizeof(int));
      ervArray[idx_integrated_variants]->alleleCnt = tmp_cnt_integrated_allele;

      int tmp_idx_allele = 0;
      for (int i = 0; i < rv_alleleCnt(rv_tmp); i++) {
        if (ifCanIntegrateAllele(rv_tmp, i, pos_start, pos_end) == true) {
          ervArray[idx_integrated_variants]->alleleIdx[tmp_idx_allele++] = i;
        }
      }
      cnt_integrated_alleles += tmp_cnt_integrated_allele;
      idx_integrated_variants++;
    } else {
      if (rv_pos(rv_tmp) > pos_end) {
        break;
      }
    }
    rv_tmp = next_RecVcf_bplus(rv_tmp);
  }

  // ------------------- Save integrated kmers into hash table -----------------
  // Calculate combinations of variants
  int *ervIdxes = (int *)calloc(cnt_integrated_variants, sizeof(int));
  for (int i = 0; i < cnt_integrated_variants; i++) ervIdxes[i] = i;
  for (int i = 1; i < cnt_integrated_variants + 1; i++) {
    Combinations *cbs =
        calculate_combinations(ervIdxes, cnt_integrated_variants, i);
    // Calculate combinations of alleles in a combination of variants
    for (int j = 0; j < cbs->cnt; j++) {
      Combinations_alleles *acbs = calculate_combinations_alleles(
          ervArray, cnt_integrated_alleles, cbs->combis[j], i);

      // Iterate all combinations of alleles in a combination of variants
      for (int k = 0; k < acbs->cnt; k++) {
        generateKmers_integrated(
            ervArray, acbs->combi_rv, acbs->combis_allele[k], acbs->length,
            lbound, rbound, pos_start, pos_end, hashTable, chrom, gv, fp_op);
      }

      for (int j = 0; j < acbs->cnt; j++) {
        free(acbs->combis_allele[j]);
      }
      free(acbs->combis_allele);
      destroy_combinations_alleles(acbs);
    }
    for (int j = 0; j < cbs->cnt; j++) {
      free(cbs->combis[j]);
    }
    free(cbs->combis);
    destroy_combinations(cbs);
  }

  // statistics_kmerHashTable(hashTable);

  // ----------------------- Write kmers into output file ----------------------
  Iterator_Kmerhash *it = init_iterator_Kmerhash(hashTable);

  while (iterator_Kmerhash_next(it) != false) {
    fprintf(fp_op, "[%" PRId32 ",%" PRId32 ",%s]\n", id_chrom,
            iterator_Kmerhash_pos(it), iterator_Kermhash_string(it));
  }

  destroy_iterator_Kmerhash(it);

  // ----------------------------- free structures -----------------------------
  destroy_kmerHashTable(hashTable);
  for (int i = 0; i < cnt_integrated_variants; i++) {
    free(ervArray[i]->alleleIdx);
    free(ervArray[i]);
  }
  free(ervArray);
  free(ervIdxes);

  return;
}

void check_files_generateKmers(Options *opts) {
  if (getAuxFile(opts) == NULL) {
    fprintf(stderr,
            "Error: argument incomplete - lack input auxiliary file.\n");
    exit(EXIT_FAILURE);
  }
  if (getVcfFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete - lack input files *.vcf/bcf\n");
    exit(EXIT_FAILURE);
  }
  if (getFaFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete - lack input files *.fa/fna/fasta\n");
    exit(EXIT_FAILURE);
  }
  if (getOutputFile(opts) == NULL) {
    fprintf(stderr,
            "Warning: output file not set. Set as default output file: %s\n",
            default_outputFile);
    setOutputFile(opts, default_outputFile);
  }
  if (opt_get_kmerLength(opts) <= 2) {
    fprintf(stderr,
            "Warning: length for generated kmer not set or invalid. Set as "
            "default: %d\n",
            default_kmerLength);
    opt_set_kmerLength(default_kmerLength, opts);
  }
}

void generateKmers(Options *opts) {
  check_files_generateKmers(opts);

  kmerLength = opt_get_kmerLength(opts);

  clock_t time_start = 0;
  clock_t time_end = 0;
  // Init structures (data storage and access)
  time_start = clock();
  GenomeFa *gf = genomeFa_loadFile(getFaFile(opts));
  time_end = clock();
  printf("... %s loaded. time: %fs\n", getFaFile(opts),
         time_convert_clock2second(time_start, time_end));
  time_start = clock();
  GenomeVcf_bplus *gv = genomeVcf_bplus_loadFile(getVcfFile(opts), 7, 6);
  time_end = clock();
  printf("... %s loaded. time: %fs\n", getVcfFile(opts),
         time_convert_clock2second(time_start, time_end));

  // Open output file
  FILE *fp_op = fopen(getOutputFile(opts), "w");
  if (fp_op == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", getOutputFile(opts));
    exit(EXIT_FAILURE);
  }

  FILE *fp_aux = fopen(getAuxFile(opts), "r");
  if (fp_aux == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", getAuxFile(opts));
    exit(EXIT_FAILURE);
  }

  uint32_t id_chrom = 0;
  uint32_t pos_start = 0;
  uint32_t pos_end = 0;
  while (fscanf(fp_aux, "[%" PRIu32 ",%" PRIu32 ",%" PRIu32 "]\n", &id_chrom,
                &pos_start, &pos_end) == 3) {
    printf("aux record: [%" PRIu32 ",%" PRIu32 ",%" PRIu32 "]\n", id_chrom,
           pos_start, pos_end);
    fprintf(fp_op, "# [%" PRIu32 ",%" PRIu32 ",%" PRIu32 "]\n", id_chrom,
            pos_start, pos_end);
    generateKmers_process(id_chrom, pos_start, pos_end, gf, gv, fp_op);
  }

  fclose(fp_op);
  fclose(fp_aux);

  // Free structures

  destroy_GenomeFa(gf);
  destroy_GenomeVcf_bplus(gv);

  return;
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

void _testSet_generateKmers() {
  // usage_auxFile();
  return;
}

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

  // Expanded interval (also ignore the kmers at the ends of the chrom)
  int32_t lbound = pos_start - (kmerLength - 1);
  lbound = lbound > 1 ? lbound : 2;
  int32_t rbound = pos_end + (kmerLength - 1);
  rbound = rbound < length_chrom ? rbound : length_chrom - 1;

  // ------------------------ Create hash table for kmers ----------------------
  // Maybe a larger size of the table would be better, but this one can do
  int32_t tableSize = (rbound - lbound) * 2;
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
      kmerHashTable_add(kmer,
                        genomeFa_absolutePos(id_chrom, lbound + tmp_offset, gf),
                        inputChar, outputChar, hashTable);
      // printf("kmer: %s, pos: %d\n", kmer, (int)(pos_start + tmp_offset));
    }
    tmp_offset++;
  }
  free(seq_ref);

  // ------------- Find all inter-variant kmers within the interval ------------

  RecVcf_bplus *rv_tmp =
      genomeVcf_bplus_getRecAfterPos(gv, name_chrom, pos_start);
  while (rv_tmp != NULL) {
    int cnt_integrated_allele = 0;
    // Check the varaint and find its inter-variant kmers
    for (int i = 0; i < rv_alleleCnt(rv_tmp); i++) {
      if (ifCanIntegrateAllele(rv_tmp, i, pos_start, pos_end) == true) {
        // genomeVcf_bplus_printRec(gv, rv_tmp);
        int32_t pos_var = rv_pos(rv_tmp);
        const char *allele_ref = rv_allele(rv_tmp, 0);
        const char *allele_alt = rv_allele(rv_tmp, i);
        const int length_ref = strlen(allele_ref);
        const int length_alt = strlen(allele_alt);

        // Divide the (l, m, r) parts and find all kmers within the interval
        // Extract (k+2) mer instead of kmer
        int32_t pos_kmer = pos_var - (kmerLength + 2 - 1);
        pos_kmer = pos_kmer >= 1 ? pos_kmer : 1;
        int32_t local_pos_alt = 0;  // 0-based pos in ALT
        while (pos_kmer != pos_var || local_pos_alt < length_alt) {
          char inputChar = 'N';
          char outputChar = 'N';
          char buf_kmer[kmerLength + 3];
          memset(buf_kmer, 0, kmerLength + 3);
          int idx_buf = 0;
          // Construct a single kmer that has an lpart
          // ----------------- lpart
          char *lpart = getSeqFromChromFa(pos_kmer, pos_var - 1, chrom);
          idx_buf = idx_buf + pos_var - pos_kmer;
          if (lpart != NULL) {  // lpart exists
            strcat(buf_kmer, lpart);
          }
          // ----------------- midpart
          int32_t pos_ref = 0;  // position of base on the reference sequence
          int idx_alt = 0;
          // Extract bases within ALT
          while (idx_buf < kmerLength + 2 &&
                 idx_alt + local_pos_alt < length_alt) {
            buf_kmer[idx_buf++] = allele_alt[idx_alt + local_pos_alt];
            idx_alt++;
          }
          // Pass bases on the reference
          pos_ref = pos_var + length_ref;  // position of next base after REF

          // ----------------- rpart
          // Extract bases after REF
          while (idx_buf < kmerLength + 2 && pos_ref <= length_chrom) {
            buf_kmer[idx_buf++] = charOfBase(getBase(chrom, pos_ref++));
          }

          // Add kmer into hash table
          // printf("kmer: %s, input char: %c, output char: %c, pos: %d\n",
          //        buf_kmer, inputChar, outputChar, pos_kmer);
          // Extract original kmer from (k+2) mer.
          char *kmer = subStr(buf_kmer, 1, (kmerLength + 2) - 1 - 1);
          inputChar = buf_kmer[0];
          outputChar = buf_kmer[kmerLength + 2 - 1];
          int32_t pos_abs = genomeFa_absolutePos(id_chrom, pos_kmer + 1, gf);
          // printf("fixed - kmer: %s, input char: %c, output char: %c, pos:
          // %d\n",
          //        kmer, inputChar, outputChar, pos_abs);
          kmerHashTable_add(kmer, pos_abs, inputChar, outputChar, hashTable);

          free(kmer);

          if (pos_kmer == pos_var) {
            local_pos_alt++;
          }
          if (lpart != NULL) {
            pos_kmer++;
          }

          free(lpart);
        }
        cnt_integrated_allele++;
      }
    }  // end of for(;;)

    // Continue iteration of variants
    if (cnt_integrated_allele == 0) {
      break;
    } else {
      if (rv_pos(rv_tmp) > pos_end) {
        break;
      }
    }
    rv_tmp = next_RecVcf_bplus(rv_tmp);
  }

  // statistics_kmerHashTable(hashTable);

  // ----------------------- Write kmers into output file ----------------------
  Iterator_Kmerhash *it = init_iterator_Kmerhash(hashTable);

  while (iterator_Kmerhash_next(it) != false) {
    fprintf(fp_op, "[%" PRId32 " %" PRId32 " %s %c %c ]\n", id_chrom,
            iterator_Kmerhash_pos(it), iterator_Kmerhash_string(it),
            iterator_KmerHash_inputChar(it), iterator_KmerHash_outputChar(it));
  }

  destroy_iterator_Kmerhash(it);

  // ----------------------------- free structures -----------------------------
  destroy_kmerHashTable(hashTable);

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
    printf("processing aux record: [%" PRIu32 ",%" PRIu32 ",%" PRIu32 "]\n",
           id_chrom, pos_start, pos_end);
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

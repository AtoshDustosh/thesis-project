#include "grbvOperations.h"

/*********************************************************************
 *                Auxiliary Structures and Functions
 ********************************************************************/

#define VAR_SNP 1
#define VAR_SMALL_INS 2
#define VAR_SMALL_DEL 3
#define VAR_MNP 4
#define VAR_SV_INS 5
#define VAR_SV_DEL 6
#define VAR_tagged 7

static int32_t sv_min_len = 0;

typedef struct VariantStatistics VariantStatistics;

/**
 * @brief  Singly-linked-list without empty header.
 */
struct VariantStatistics {
  int32_t type_var;
  char *tag;  // NULL for variants that are not tagged.
  int32_t cnt;
  VariantStatistics *next;
};

/**
 * @brief  Initialize an object of VariantStatistics.
 */
VariantStatistics *init_VariantStatistics(int32_t type_var, const char *tag) {
  VariantStatistics *vs =
      (VariantStatistics *)malloc(sizeof(VariantStatistics));
  vs->type_var = type_var;
  if (type_var == VAR_tagged) {
    vs->tag = strdup(tag);
  } else {
    vs->tag = NULL;
  }
  vs->cnt = 0;
  vs->next = NULL;
  return vs;
}

void print_VaraintStatistics(VariantStatistics *vs) {
  char *string_type_var = NULL;
  switch (vs->type_var) {
    case VAR_SNP: {
      string_type_var = "snp";
      break;
    }
    case VAR_SMALL_INS: {
      string_type_var = "small_ins";
      break;
    }
    case VAR_SMALL_DEL: {
      string_type_var = "small_del";
      break;
    }
    case VAR_MNP: {
      string_type_var = "mnp";
      break;
    }
    case VAR_SV_INS: {
      string_type_var = "sv_ins";
      break;
    }
    case VAR_SV_DEL: {
      string_type_var = "sv_del";
      break;
    }
    case VAR_tagged: {
      string_type_var = vs->tag;
      break;
    }
    default: {
      fprintf(stderr,
              "Error: unexpected error when printing statistics for a kind of "
              "variant. \n");
      exit(EXIT_FAILURE);
    }
  }
  printf("(%s): %" PRId32 "\n", string_type_var, vs->cnt);
}

/**
 * @brief  Destroy an object of VariantStatistics and return the next object
 * linked to it.
 */
VariantStatistics *destroy_VariantStatistics(VariantStatistics *vs) {
  VariantStatistics *next = vs->next;
  free(vs->tag);
  free(vs);
  return next;
}

/**
 * @brief  Update statistics. If there is no such type of variant yet, do
 * nothing. Otherwise, update the statistics for the variant's type.
 */
static void update_VariantStatistics(int32_t type_var, const char *tag,
                                     VariantStatistics *vs) {
  VariantStatistics *tmp_vs = NULL;
  while (vs != NULL) {
    if (vs->type_var == type_var) {
      if (vs->type_var != VAR_tagged) {  // same type of variant
        vs->cnt++;
        return;
      } else {                            // if this is a tagged type
        if (strcmp(tag, vs->tag) == 0) {  // same tag
          vs->cnt++;
          return;
        } else {  // different tags
          tmp_vs = vs;
          vs = vs->next;
        }
      }
    } else {
      tmp_vs = vs;
      vs = vs->next;
    }
  }
  VariantStatistics *new_vs = init_VariantStatistics(type_var, tag);
  tmp_vs->next = new_vs;
  new_vs->cnt++;
}

void process_VariantStatistics(bcf1_t *rec, VariantStatistics *vs) {
  assert(rec->n_allele >= 1);
  const char *ref = rec->d.allele[0];
  const int length_ref = strlen(ref);
  for (int i = 1; i < rec->n_allele; i++) {
    int32_t type_var = 0;
    const char *tag;
    const char *alt = rec->d.allele[i];
    const int length_alt = strlen(alt);

    if (alt[0] == '<') {  // Tagged
      type_var = VAR_tagged;
      tag = alt;
    } else {
      if (length_ref == 1) {
        if (length_alt == 1) {  // SNP
          type_var = VAR_SNP;
        } else {
          if (length_alt >= sv_min_len) {  // sv_ins
            type_var = VAR_SV_INS;
          } else {  // small_ins
            type_var = VAR_SMALL_INS;
          }
        }
      } else {
        if (length_alt == 1) {
          if (length_ref >= sv_min_len) {  // sv_del
            type_var = VAR_SV_DEL;
          } else {  // small_del
            type_var = VAR_SMALL_DEL;
          }
        } else {
          type_var = VAR_MNP;
        }
      }
    }
    update_VariantStatistics(type_var, tag, vs);
  }
  return;
}

/*********************************************************************
 *                           GRBV operations
 ********************************************************************/

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

void statistics_vcf(Options *opts) {
  // Check input parameters
  sv_min_len = getSVminLen(opts);
  const char *filePath_vcf = getVcfFile(opts);
  htsFile *fp_in = NULL;
  if (filePath_vcf == NULL) {
    fprintf(stderr,
            "Error: vcf file unspecified yet. Cannot collect statistics. \n");
    exit(EXIT_FAILURE);
  } else {
    fp_in = hts_open(filePath_vcf, "r");
  }
  const char *filePath_out = getOutputFile(opts);
  FILE *fp_out = NULL;
  if (filePath_out == NULL) {
    fprintf(stderr,
            "Warning: output file unspecified. Set as default. Output will be "
            "printed to the console.\n");
    fp_out = stdout;
  } else {
    fp_out = fopen(filePath_out, "w");
  }

  // Collect statistics
  bcf_hdr_t *hdr = bcf_hdr_read(fp_in);
  bcf1_t *rec = bcf_init1();

  if (hdr == NULL) {
    fprintf(stderr, "Error: failed creating bcf header struct.\n");
    exit(EXIT_FAILURE);
  }
  if (rec == NULL) {
    fprintf(stderr, "Error: memory not enough for new bcf1_t object.\n");
    exit(EXIT_FAILURE);
  }

  VariantStatistics *vs = init_VariantStatistics(VAR_SNP, NULL);

  while (bcf_read1(fp_in, hdr, rec) >= 0) {
    // The "bcf_unpack" method must be called for every new bcf1_t object
    bcf_unpack(rec, BCF_UN_ALL);
    bcf1_t *tmpRec = bcf_dup(rec);
    bcf_unpack(tmpRec, BCF_UN_ALL);
    process_VariantStatistics(tmpRec, vs);
  }

  while (vs != NULL) {
    char *string_type_var = NULL;
    switch (vs->type_var) {
      case VAR_SNP: {
        string_type_var = "snp";
        break;
      }
      case VAR_SMALL_INS: {
        string_type_var = "small_ins";
        break;
      }
      case VAR_SMALL_DEL: {
        string_type_var = "small_del";
        break;
      }
      case VAR_MNP: {
        string_type_var = "mnp";
        break;
      }
      case VAR_SV_INS: {
        string_type_var = "sv_ins";
        break;
      }
      case VAR_SV_DEL: {
        string_type_var = "sv_del";
        break;
      }
      case VAR_tagged: {
        string_type_var = vs->tag;
        break;
      }
      default: {
        fprintf(
            stderr,
            "Error: unexpected error when printing statistics for a kind of "
            "variant. \n");
        exit(EXIT_FAILURE);
      }
    }
    fprintf(fp_out, "(%s): %" PRId32 "\n", string_type_var, vs->cnt);
    VariantStatistics *next = destroy_VariantStatistics(vs);
    vs = next;
  }

  // Close file streams
  if (fp_in != NULL) {
    hts_close(fp_in);
  }
  if (fp_out != stdout) {
    fclose(fp_out);
  }
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

void _testSet_grbvOperations() { return; }
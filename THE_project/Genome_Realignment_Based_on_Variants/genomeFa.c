#include "genomeFa.h"

/*********************************************************************
 *                Auxiliary Structures and Functions
 ********************************************************************/

/**
 * @brief  A linked list without empty header. Designed for initializing
 * memories before loading reference genome.
 */
typedef struct ProfileFa ProfileFa;

struct ProfileFa {
  char *info;                // info field of the chromosome
  uint32_t length;           // length of the chromosome / cnt of bases
  uint32_t absolute_offset;  // offset for absolute position of this chromosome
  ProfileFa *next;           // next chromosome's profile, if exists
};

// *************************** Declarations **************************

/**
 * @brief  Initialize profiles from an input reference genome file and return
 * pointer to a linked list of ProfileFas without empty header.
 */
ProfileFa *init_ProfileFa(const char *filePath);

void print_ProfileFa(ProfileFa *profile);

/**
 * @brief  Destroy an object of ProfileFa and return pointer to the profile of
 * next chromosome.
 */
ProfileFa *destroy_Profile(ProfileFa *profile);

// ************************ Implementations **************************

ProfileFa *init_ProfileFa(const char *filePath) {
  FILE *fp = fopen(filePath, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: empty file pointer. \n");
    exit(EXIT_FAILURE);
  }

  ProfileFa *profile = (ProfileFa *)calloc(1, sizeof(ProfileFa));
  ProfileFa *tmp_profile = profile;  // the last created profile
  ProfileFa *new_profile = NULL;

  uint32_t pos_bp = 1;  // 1-based position within chromosome

  char buf_info[MAX_INFO_LENGTH];  // buffer for info line of the file
  int idx_buf_info = 0;            // index for buf_info

  bool isBpLine = false;

  char tmpCh;
  while ((tmpCh = fgetc(fp)) != EOF) {
    if (tmpCh == '>') {
      if (pos_bp != 1) {
        new_profile->length = pos_bp - 1;
        new_profile->absolute_offset =
            tmp_profile->absolute_offset + tmp_profile->length;
        tmp_profile->next = new_profile;
        tmp_profile = new_profile;
      }
      isBpLine = false;
      idx_buf_info = 0;
      pos_bp = 1;
      buf_info[idx_buf_info++] = tmpCh;
      continue;
    }
    // Extract the profile
    if (isBpLine == false) {  // Extract info line
      if (tmpCh == '\n') {    // Reach the end of info line
        buf_info[idx_buf_info] = '\0';
        new_profile = (ProfileFa *)calloc(1, sizeof(ProfileFa));
        new_profile->info = strdup(buf_info);
        isBpLine = true;
      } else {  // Buffer character received from info line
        buf_info[idx_buf_info++] = tmpCh;
      }
    } else {  // Extract bases
      if (tmpCh == '\n') {
        continue;
      } else {
        pos_bp++;
      }
    }
  }
  if (pos_bp != 1) {
    new_profile->length = pos_bp - 1;
    new_profile->absolute_offset =
        tmp_profile->absolute_offset + tmp_profile->length;
    tmp_profile->next = new_profile;
    tmp_profile = new_profile;
  }

  fclose(fp);

  return destroy_Profile(profile);
}

void print_ProfileFa(ProfileFa *profile) {
  ProfileFa *tmp_profile = profile;
  while (tmp_profile != NULL) {
    printf("info: %s, length: %" PRIu32 "\n", tmp_profile->info,
           tmp_profile->length);
    tmp_profile = tmp_profile->next;
  }
}

ProfileFa *destroy_Profile(ProfileFa *profile) {
  ProfileFa *next = profile->next;
  free(profile->info);
  free(profile);
  return next;
}

/**
 * @brief  (helper function) Code the bpBuf array (a string with A,C,G,T) into a
 * uint64_t integer. Should not be called by the user. The calculation's
 * purpose: "ACGT" -> (0b) 001 010 011 100 000 000 ...... 000 0. Note that bases
 * start from the left and there may be some "0"s left out not filled. The bpBuf
 * string must contain no more than BP_PER_UINT64 chars (bases).
 */
static uint64_t codeBpBuf(char *bpBuf) {
  uint64_t codedBp = 0x0;
  int i = 0;
  while (i != BP_PER_UINT64 && *bpBuf != '\0') {
    /*
     * The following calculation's purpose:
     * "ACGT" -> (0b) 001 010 011 100 000 000 ...... 000 0
     * The bpBuf string must contain no more than BP_PER_UINT64 chars (bases).
     */
    codedBp =
        codedBp | (((uint64_t)baseOfChar(*bpBuf))
                   << (sizeof(uint64_t) * 8) - BASE_CODE_LENGTH * (i + 1));
    bpBuf++;
    i++;
  }
  return codedBp;
}

/*********************************************************************
 *                     Structure Declarations
 ********************************************************************/

// This is actually a linked-list with an empty header
struct ChromFa {
  uint64_t *codedBases;      // binary coded bases using uint64_t array
  uint32_t length;           // length of chrom / number of bases (uncoded)
  uint32_t absolute_offset;  // offset for absolute position of this chromosome
  char *info;                // info of chrom
  char *name;
  struct ChromFa *next;
};

struct GenomeFa {
  uint16_t chromCnt;
  ChromFa *chroms;
};

/**
 * @brief Initialize a ChromFa object and return the pointer to it. The
 * successfully returned object must be destroyed later using destroy_ChromFa()
 */
static ChromFa *init_ChromFa() {
  ChromFa *cf = (ChromFa *)calloc(1, sizeof(ChromFa));
  if (cf == NULL) {
    fprintf(stderr, "Error: no enough memory for a new chrom.\n");
    exit(EXIT_FAILURE);
  }
  return cf;
}

static ChromFa *init_ChromFa_info_length(const char *info,
                                         const uint32_t length) {
  ChromFa *cf = (ChromFa *)calloc(1, sizeof(ChromFa));
  if (cf == NULL) {
    fprintf(stderr, "Error: no enough memory for a new chrom.\n");
    exit(EXIT_FAILURE);
  }
  uint32_t length_info = strlen(info);
  cf->info = strdup(info);
  cf->length = length;
  cf->codedBases =
      (uint64_t *)calloc((length - 1) / BP_PER_UINT64 + 1, sizeof(uint64_t));

  // Find "name" of the chromosome (">name other")
  if (info[0] != '>') {
    fprintf(stderr, "Error: info line of chromosome invalid - %s\n", info);
    exit(EXIT_FAILURE);
  }
  char buf_name[MAX_INFO_LENGTH];
  memset(buf_name, 0, MAX_INFO_LENGTH);
  for (int i = 0; i < length_info; i++) {
    char tmp_char = info[i + 1];
    if (tmp_char == ' ' || tmp_char == '\t') {
      break;
    } else {
      buf_name[i] = tmp_char;
    }
  }
  cf->name = strdup(buf_name);
  return cf;
}

static ChromFa *init_ChromFa_profile(ProfileFa *profile) {
  ChromFa *cf = (ChromFa *)calloc(1, sizeof(ChromFa));
  if (cf == NULL) {
    fprintf(stderr, "Error: no enough memory for a new chrom.\n");
    exit(EXIT_FAILURE);
  }
  const char *info = profile->info;
  const int length_info = strlen(info);
  const uint32_t length = profile->length;
  const uint32_t absolute_offset = profile->absolute_offset;

  cf->info = strdup(info);
  cf->length = length;
  cf->absolute_offset = absolute_offset;
  cf->codedBases =
      (uint64_t *)calloc((length - 1) / BP_PER_UINT64 + 1, sizeof(uint64_t));

  // Find "name" of the chromosome (">name other")
  if (info[0] != '>') {
    fprintf(stderr, "Error: info line of chromosome invalid - %s\n", info);
    exit(EXIT_FAILURE);
  }
  char buf_name[MAX_INFO_LENGTH];
  memset(buf_name, 0, MAX_INFO_LENGTH);
  for (int i = 0; i < length_info; i++) {
    char tmp_char = info[i + 1];
    if (tmp_char == ' ' || tmp_char == '\t') {
      break;
    } else {
      buf_name[i] = tmp_char;
    }
  }
  cf->name = strdup(buf_name);
  return cf;
}

/**
 * @brief Destroy a ChromFa object. Should not be called by the user.
 */
static void destroy_ChromFa(ChromFa *cf) {
  if (cf == NULL) {
    fprintf(stderr,
            "Error: null pointer occurred when destroying a ChromFa object\n");
    exit(EXIT_FAILURE);
  }
  free(cf->codedBases);
  free(cf->info);
  free(cf->name);
  free(cf);
}

/**
 * @brief Add a ChromFa object into the GenomeFa object (linked to the end of
 * the linked-list with an empty header).
 */
static void addChromToGenome(ChromFa *cf, GenomeFa *gf) {
  if (cf == NULL || gf == NULL) {
    fprintf(stderr, "Error: null pointer occurred for ChromFa or GenomeFa\n");
    exit(EXIT_FAILURE);
  }
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    if (tmpCf->next == NULL) {
      tmpCf->next = cf;
      gf->chromCnt++;
      return;
    }
    tmpCf = tmpCf->next;
  }
}

GenomeFa *init_GenomeFa() {
  GenomeFa *gf = (GenomeFa *)malloc(sizeof(GenomeFa));
  if (gf == NULL) {
    fprintf(stderr, "Error: no enough memory for a new GenomeFa object.\n");
    exit(EXIT_FAILURE);
  }
  gf->chromCnt = 0;
  gf->chroms = init_ChromFa();
  gf->chroms->codedBases = (uint64_t *)calloc(1, sizeof(uint64_t));
  gf->chroms->info = (char *)calloc(strlen(">header chrom") + 1, sizeof(char));
  gf->chroms->name = (char *)calloc(strlen("header") + 1, sizeof(char));
  strcpy(gf->chroms->info, ">header chrom");
  strcpy(gf->chroms->name, "header");
  return gf;
}

void destroy_GenomeFa(GenomeFa *gf) {
  if (gf == NULL) {
    fprintf(stderr,
            "Error: null pointer occurred when destroying a GenomeFa object\n");
    exit(EXIT_FAILURE);
  }
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    ChromFa *nxtCf = tmpCf->next;
    destroy_ChromFa(tmpCf);
    tmpCf = nxtCf;
  }
  free(gf);
}

/*********************************************************************
 *                           Data Accessor
 ********************************************************************/

inline uint32_t chromFa_length(ChromFa *cf) { return cf->length; }

inline const char *chromFa_name(ChromFa *cf) { return cf->name; }

ChromFa *getChromFromGenomeFabyInfo(const char *info, GenomeFa *gf) {
  if (gf == NULL) {
    fprintf(stderr, "Error: null pointer occurred for GenomeFa\n");
    assert(false);
  }
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    char *tmpInfo = tmpCf->info;
    if (tmpInfo != NULL) {
      if (strcmp(tmpInfo, info) == 0) {
        return tmpCf;
      }
    }
    tmpCf = tmpCf->next;
  }
  // return NULL if didn't found any chrom with the same info field
  return NULL;
}

ChromFa *getChromFromGenomeFabyName(const char *name, GenomeFa *gf) {
  if (gf == NULL || name == NULL) {
    fprintf(stderr, "Error: null pointer occurred for GenomeFa\n");
    assert(false);
  }
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    char *tmpName = tmpCf->name;
    if (tmpName != NULL) {
      if (strcmp(tmpName, name) == 0) {
        return tmpCf;
      }
    }
    tmpCf = tmpCf->next;
  }
  return NULL;
}

ChromFa *getChromFromGenomeFabyIndex(uint32_t idx, GenomeFa *gf) {
  if (gf == NULL) {
    fprintf(stderr, "Error: null pointer occurred for GenomeFa\n");
    assert(false);
  }
  ChromFa *tmpCf = gf->chroms;
  for (uint32_t i = 0; i < idx; i++) {
    tmpCf = tmpCf->next;
  }

  return tmpCf;
}

Base getBase(ChromFa *cf, uint32_t pos) {
  // Spent some time for this. Bug occurred. Bug located. Bug fixed ...
  uint32_t arrayLength = (cf->length - 1) / BP_PER_UINT64 + 1;
  uint32_t arrayIdx = (pos - 1) / BP_PER_UINT64;

  if (arrayIdx >= arrayLength || pos <= 0) {
    return BASE_INVALID;
  } else {
    // The following commented lines were used for debugging.
    uint64_t arrayElement = cf->codedBases[arrayIdx];
    // 1 % 21 = 1, 21 % 21 = 0. This is not good. So (21 - 1) % 21 + 1 = 21.
    // (pos, baseIdx) examples: (1,1), (21,21), (22,1), (43,1)
    uint8_t baseIdx = (pos - 1) % BP_PER_UINT64 + 1;
    // (baseIdx, offset) examples, (21,1), (20,4), (1, 61)
    uint8_t offset = 1 + BASE_CODE_LENGTH * (BP_PER_UINT64 - baseIdx);

    Base retVal = (Base)((arrayElement >> offset) & BASE_MASK_RIGHT);
    return retVal;
    // return (Base)(
    //     (cf->codedBases[arrayIdx] >>
    //          (1 + BASE_CODE_LENGTH *
    //                   (BP_PER_UINT64 - ((pos - 1) % BP_PER_UINT64 + 1))) &
    //      BASE_MASK_RIGHT));
  }
}

char *getSeqFromChromFa(int64_t start, int64_t end, ChromFa *cf) {
  if (cf == NULL || end < start) {
    // assert(fprintf(stderr, "Report: invalid parameters for getSeq.\n") > 0);
    return NULL;
  }
  uint64_t seqLength = (end - start + 1);
  //  "+1" is in consideration of '\0' at the end of a string.
  char *seq = (char *)calloc(seqLength + 1, sizeof(char));

  for (int i = 0; i < seqLength; i++) {
    seq[i] = charOfBase(getBase(cf, start + i));
  }
  seq[seqLength] = '\0';

  return seq;
}

int32_t genomeFa_absolutePos(uint32_t id_chrom, int32_t pos, GenomeFa *gf) {
  if (gf == NULL) {
    fprintf(stderr, "Error: null pointer occurred for GenomeFa\n");
    assert(false);
  }
  ChromFa *tmpCf = gf->chroms;
  for (uint32_t i = 0; i < id_chrom; i++) {
    tmpCf = tmpCf->next;
  }
  return pos + tmpCf->absolute_offset;
}

/*********************************************************************
 *                      Data Loading and Writing
 ********************************************************************/

GenomeFa *genomeFa_loadFile(char *filePath) {
  GenomeFa *gf = init_GenomeFa();
  // Initialie profiles of all chromosomes

  // Add chromosomes according to profiles of them and destroy profiles
  ProfileFa *tmp_profile = init_ProfileFa(filePath);
  while (tmp_profile != NULL) {
    ChromFa *tmp_cf = init_ChromFa_profile(tmp_profile);
    addChromToGenome(tmp_cf, gf);
    // printf("absolute offsets: %" PRIu32 "\n", tmp_profile->absolute_offset);
    tmp_profile = destroy_Profile(tmp_profile);
  }

  char buf_bases[BP_PER_UINT64 + 1];
  int idx_buf_bases = 0;

  bool isBpLine = false;
  uint32_t pos_bp = 1;  // 1-based position of bases within chrom

  ChromFa *current_cf = gf->chroms;  // Temporary under processing chrom

  FILE *fp = fopen(filePath, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open file - %s\n", filePath);
    exit(EXIT_FAILURE);
  }

  char tmp_char;
  while ((tmp_char = fgetc(fp)) != EOF) {
    if (tmp_char == '>') {
      if (idx_buf_bases != 0) {
        buf_bases[idx_buf_bases] = '\0';
        current_cf->codedBases[(pos_bp - 1) / BP_PER_UINT64] =
            codeBpBuf(buf_bases);
        idx_buf_bases = 0;
      }
      isBpLine = false;
      continue;
    }
    if (isBpLine == false) {
      if (tmp_char == '\n') {  // reach the end of info line
        current_cf = current_cf->next;
        isBpLine = true;
        pos_bp = 1;
      } else {  // Ignore the info line characters
        continue;
      }
    } else {  // Code bases into 64-bit integers and store them into chromsomes
      if (tmp_char == '\n') {  // Ignore the line shifters
        continue;
      } else {
        buf_bases[idx_buf_bases++] = tmp_char;
        if (idx_buf_bases == BP_PER_UINT64) {  // when the buf_bases is full
          buf_bases[idx_buf_bases] = '\0';
          current_cf->codedBases[(pos_bp - 1) / BP_PER_UINT64] =
              codeBpBuf(buf_bases);
          idx_buf_bases = 0;
        }
        pos_bp++;
      }
    }
  }
  // If the file doesn't have an empty line at the end, flush the buf_bases
  if (idx_buf_bases != 0) {
    buf_bases[idx_buf_bases] = '\0';
    current_cf->codedBases[(pos_bp - 1) / BP_PER_UINT64] = codeBpBuf(buf_bases);
    idx_buf_bases = 0;
  }

  fclose(fp);

  // printGenomeFa(gf);
  return gf;
}

void writeGenomeFaIntoFile(GenomeFa *gf, const char *filePath) {
  FILE *fp = fopen(filePath, "w");
  ChromFa *tmpCf = gf->chroms->next;

  while (tmpCf != NULL) {
    fprintf(fp, "%s\n", tmpCf->info);
    uint32_t baseCnt = tmpCf->length;
    int i = 1;
    for (i = 1; i <= baseCnt; i++) {
      char tmpBp = charOfBase(getBase(tmpCf, i));
      fprintf(fp, "%c", tmpBp);
      if (i % BP_PER_LINE == 0) {
        fprintf(fp, "\n");
      }
    }
    if (i % BP_PER_LINE != 0) {
      fprintf(fp, "\n");
    }
    tmpCf = tmpCf->next;
  }
  fclose(fp);
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

static void _test_StructureLinks() {
  /*
   * Please note that this test is not the valid usage of the APIs.
   * Normally you should not access the fields of the structures directly.
   */
  GenomeFa *gf = init_GenomeFa();
  ChromFa *cf1 = init_ChromFa();
  ChromFa *cf2 = init_ChromFa();
  ChromFa *cf3 = init_ChromFa();
  ChromFa *cf4 = init_ChromFa();
  ChromFa *cf5 = init_ChromFa();

  cf1->length = 1;
  cf1->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 1);
  cf1->codedBases[0] = 0x92449244;
  cf1->info = (char *)malloc(sizeof(char *) * strlen(">chr1"));
  strcpy(cf1->info, ">chr1");
  addChromToGenome(cf1, gf);

  cf2->length = 2;
  cf2->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 2);
  cf2->codedBases[0] = 0x92449244;
  cf2->codedBases[1] = 0x92449244;
  cf1->info = (char *)malloc(sizeof(char *) * strlen(">chr2"));
  strcpy(cf1->info, ">chr2");
  addChromToGenome(cf2, gf);

  cf3->length = 3;
  cf3->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 3);
  cf3->codedBases[0] = 0x92449244;
  cf3->codedBases[1] = 0x92449244;
  cf3->codedBases[2] = 0x92449244;
  cf1->info = (char *)malloc(sizeof(char *) * strlen(">chr3"));
  strcpy(cf1->info, ">chr3");
  addChromToGenome(cf3, gf);

  cf4->length = 4;
  cf4->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 4);
  cf4->codedBases[0] = 0x92449244;
  cf4->codedBases[1] = 0x92449244;
  cf4->codedBases[2] = 0x92449244;
  cf4->codedBases[3] = 0x92449244;
  cf1->info = (char *)malloc(sizeof(char *) * strlen(">chr4"));
  strcpy(cf1->info, ">chr4");
  addChromToGenome(cf4, gf);

  cf5->length = 5;
  cf5->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 5);
  cf5->codedBases[0] = 0x92449244;
  cf5->codedBases[1] = 0x92449244;
  cf5->codedBases[2] = 0x92449244;
  cf5->codedBases[3] = 0x92449244;
  cf5->codedBases[4] = 0x92449244;
  cf1->info = (char *)malloc(sizeof(char *) * strlen(">chr5"));
  strcpy(cf1->info, ">chr5");
  addChromToGenome(cf5, gf);

  destroy_GenomeFa(gf);
}

static void _test_CodingBases() {
  char *bases[] = {"AAAACCCCGGGGTTTTAAAAC", "TACCCCAGCAAGACCAGACAT", "ACGT"};
  uint64_t codedBases[] = {0x2494926DB9242494, 0x85248b44b2916518,
                           0x29C0000000000000};
  ChromFa *cf = init_ChromFa();
  cf->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 3);
  cf->length = 46;
  cf->info = (char *)malloc(sizeof(char *) * strlen(">codingBasesTestCf"));
  strcpy(cf->info, ">codingBasesTestCf");

  // test the calculation of coded bases
  for (int i = 0; i < (sizeof(codedBases) / sizeof(uint64_t)); i++) {
    assert(codedBases[i] == codeBpBuf(bases[i]));
    cf->codedBases[i] = codeBpBuf(bases[i]);
  }

  // test the extraction of coded bases from ChromFa
  assert((getBase(cf, 1) == BASE_A) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 1), BASE_A)) < 0);
  assert((getBase(cf, 10) == BASE_G) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 10), BASE_G)) < 0);
  assert((getBase(cf, 21) == BASE_C) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 21), BASE_C)) < 0);
  assert((getBase(cf, 22) == BASE_T) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 22), BASE_T)) < 0);
  assert((getBase(cf, 24) == BASE_C) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 24), BASE_C)) < 0);
  assert((getBase(cf, 42) == BASE_T) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 42), BASE_T)) < 0);
  assert((getBase(cf, 43) == BASE_A) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 43), BASE_A)) < 0);
  assert((getBase(cf, 44) == BASE_C) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 44), BASE_C)) < 0);
  assert((getBase(cf, 46) == BASE_T) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 46), BASE_T)) < 0);
  assert((getBase(cf, 999) == BASE_INVALID) ||
         (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                  getBase(cf, 999), BASE_INVALID)) < 0);

  destroy_ChromFa(cf);
}

static void _test_CodingBases_Plus() {
  const char *chromName[] = {"chr1", "chr2", "chr3", "chr4"};
  const char *chromSeq[] = {
      "AAAAAAAAAAAAACCCCCCTATACCCCAGCAAGACCAGACATCCCCCCCCCATAGACACCCCCACGTACGTA"
      "CGTACGTACGTCCCCCCCCGATACAAGTAGACCCCCAGCATACCATACCCCCCCCCCCCCCAAAAAAAAAAA"
      "AAAAAAAAAAAAAAA",
      "CCCCCCCCCGGGGGGGGGGGGATAACATAAGGGGATACGGGGTAGAATAGGGGGGGGGGACGTACGTACGTA"
      "CGTACGTGGGGGTAGAGAGACAGAGAGATACGGGGGGGGGGGGGGTAGAGAGGGGGGCGTACGTACGTACGT"
      "ACGTAGGGGGATACACAGATACGATTAACGATAACATATAATGGGGGGGGGGGGGCCCCCCCCCCCCCCCCC"
      "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAACCCCCCTATACC"
      "CCAGCAAGACCAGACATCCCCGGCCCATAGACACCCCCACGT",
      "GGGGGGGGGTTTTTTTTTTTATAAAATAGGCCACACCAGAAATAGACCTTTTTTTTTTTACGTACGTACGTA"
      "CGTACGTTTTTTTTTCACAGTAAGCATATTTTTTTGAGAGAGAAAATTTAATTTTTTCGTACGTACGTACGT"
      "ACGTATTTTTTTAGAATATACCCAAATATTTTGAGAATAAACCATTTTTTTTTTTGGGGGGGGGGGGGGGGG"
      "GGGGGG",
      "TTTTTTTTTAAAAAAAAACGCTAGCATCGATAGCTAGGAATATATAGCAAAAAAAAAAAACGTACGTACGTA"
      "CGTACGTAAGGGTTTAACCAAGTAACAAAAAAAAAAACCCCCGGATTTTAAAAAAAACGTACGTACGTACGT"
      "ACGTAAAAAAGGGAAATTTAACAAAAAAAAGAGAATTACCCCCAATAAAAAAAAATTTTTTTTTTTTT"};
  const uint32_t testCnt = sizeof(chromName) / sizeof(char *);
  GenomeFa *gf = genomeFa_loadFile("data/test.fa");

  for (int i = 0; i < testCnt; i++) {
    ChromFa *cf = getChromFromGenomeFabyName(chromName[i], gf);
    for (int j = 0; j < strlen(chromSeq[i]); j++) {
      char bpGet = charOfBase(getBase(cf, j + 1));
      if (bpGet != chromSeq[i][j]) {
        assert(bpGet == chromSeq[i][j]);
      }
    }
  }
  destroy_GenomeFa(gf);
}

static void _test_SeqExtractor() {
  const uint64_t testStart[] = {1, 53, 71};
  const uint64_t testEnd[] = {34, 85, 71};
  const char *testSeq[] = {"AAAAAAAAAAAAACCCCCCTATACCCCAGCAAGA",
                           "TAGACACCCCCACGTACGTACGTACGTACGTCC", "T"};
  const char *testName[] = {"chr1", "chr1", "chr1"};
  const uint32_t testCnt = sizeof(testStart) / sizeof(uint64_t);
  GenomeFa *gf = genomeFa_loadFile("data/test.fa");

  for (int i = 0; i < testCnt; i++) {
    ChromFa *tmpCf = getChromFromGenomeFabyName(testName[i], gf);
    char *seq = getSeqFromChromFa(testStart[i], testEnd[i], tmpCf);
    // printf("seq : %s\n", seq);
    // printf("test: %s\n", testSeq[i]);
    assert(strcmp(seq, testSeq[i]) == 0);
    free(seq);
  }

  destroy_GenomeFa(gf);
}

static void _test_Loader() {
  GenomeFa *gf = genomeFa_loadFile("data/example.fa");

  // printGenomeFa(gf);

  uint32_t testFaCfs_length[] = {151, 203, 152, 142};
  char *testFaCfs_info[] = {
      ">chr1  LN:151  AS:GRCh38", ">chr2  LN:203  AS:GRCh38",
      ">chr13  LN:152  AS:GRCh38", ">chr14  LN:142  AS:GRCh38"};

  ChromFa *tmpCf = gf->chroms;
  tmpCf = tmpCf->next;
  for (int i = 0; i < sizeof(testFaCfs_length) / sizeof(int); i++) {
    // printf("tmpCfLN: %d, testLN: %d\n", tmpCf->length, testFaCfs_length[i]);
    // printf("tmpCfinfo: \"%s\", testinfo: \"%s\"\n", tmpCf->info,
    //        testFaCfs_info[i]);
    assert(tmpCf->length == testFaCfs_length[i]);
    assert(strcmp(tmpCf->info, testFaCfs_info[i]) == 0);
    tmpCf = tmpCf->next;
  }

  uint64_t codedCh1[] = {0x2492492492492492, 0x2492492492492492,
                         0x2492492492492492, 0x2492492492492492,
                         0x2492492492492492, 0x2492492492492492,
                         0x2492492492492492, 0x2490000000000000};
  uint64_t codedCh2[] = {0x4924924924924924, 0x4924924924924924,
                         0x4924924924924924, 0x4924924924924924,
                         0x4924924924924924, 0x4924924924924924,
                         0x4924924924924924, 0x4924924924924924,
                         0x4924924924924924, 0x4924924924800000};
  uint64_t codedCh3[] = {0x6db6db6db6db6db6, 0x6db6db6db6db6db6,
                         0x6db6db6db6db6db6, 0x6db6db6db6db6db6,
                         0x6db6db6db6db6db6, 0x6db6db6db6db6db6,
                         0x6db6db6db6db6db6, 0x6db6000000000000};
  uint64_t codedCh4[] = {0x9249249249249248, 0x9249249249249248,
                         0x9249249249249248, 0x9249249249249248,
                         0x9249249249249248, 0x9249249249249248,
                         0x9249249249240000};
  uint64_t *codedChs[] = {codedCh1, codedCh2, codedCh3, codedCh4};
  tmpCf = gf->chroms->next;
  for (int i = 0; i < gf->chromCnt; i++) {
    for (int j = 0; j < (tmpCf->length - 1) / BP_PER_UINT64; j++) {
      assert(tmpCf->codedBases[j] == codedChs[i][j]);
    }
    tmpCf = tmpCf->next;
  }

  destroy_GenomeFa(gf);
}

static void _test_Writer() {
  GenomeFa *gf = genomeFa_loadFile("data/test.fa");

  // printGenomeFa(gf);

  writeGenomeFaIntoFile(gf, "data/output.fa");

  destroy_GenomeFa(gf);
}

static void _test_Loader_refactored() {
  GenomeFa *gf = genomeFa_loadFile("data/example.fa");

  // printGenomeFa(gf);

  destroy_GenomeFa(gf);
}

void _testSet_genomeFa() {
  _test_StructureLinks();
  _test_CodingBases();
  _test_CodingBases_Plus();
  _test_SeqExtractor();
  _test_Loader();
  _test_Writer();
  _test_Loader_refactored();
}

void printGenomeFa(GenomeFa *gf) {
  printf("number of chroms: %d\n", gf->chromCnt);
  // Ignore the empty header of the linked list
  ChromFa *tmpCf = gf->chroms->next;
  while (tmpCf != NULL) {
    printf("info: %s\n", tmpCf->info);
    printf("length: %" PRIu32 "\n", tmpCf->length);
    printf("absolute offset: %" PRIu32 "\n", tmpCf->absolute_offset);

    uint32_t baseCnt = tmpCf->length;
    int i = 0;
    for (i = 1; i <= baseCnt; i++) {
      char tmpBp = charOfBase((getBase(tmpCf, i)));
      fprintf(stdout, "%c", tmpBp);
      if (i % BP_PER_LINE == 0) {
        fprintf(stdout, "\n");
      }
    }
    if (i % BP_PER_LINE != 0) {
      fprintf(stdout, "\n");
    }
    tmpCf = tmpCf->next;
  }
}

void printChromFa(ChromFa *cf) {
  printf("info: %s\n", cf->info);
  printf("length: %" PRIu32 "\n", cf->length);

  uint32_t arrayLength = 0;
  if (cf->length != 0) {
    arrayLength = (cf->length - 1) / BP_PER_UINT64 + 1;
  }
  for (int i = 0; i < arrayLength; i++) {
    printf("\t0x%" PRIx64 "\t", cf->codedBases[i]);
    if ((i + 1) % 5 == 0) {
      printf("\n");
    }
  }
  printf("\n");
}

void printGenomeFa_brief(GenomeFa *gf) {
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    printf("%s\n", tmpCf->info);
    tmpCf = tmpCf->next;
  }
}

#include "genomeFa.h"

inline static void assert_bases_equal(Base x, Base y) {
  // note that the "||" operation can terminate the condition statement
  // instantly when "x == y" succeeds and thus printf will not print. (a great
  // way for debugging using "assert(...)")
  assert(x == y || (fprintf(stderr, "calc: 0x%" PRIx8 ", true: 0x%" PRIx8 "\n",
                            x, y) >= 0));
}

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

  strcpy(cf1->info, ">chr1");
  cf1->length = 1;
  cf1->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 1);
  cf1->codedBases[0] = 0x92449244;
  addChromToGenome(cf1, gf);

  strcpy(cf2->info, ">chr2");
  cf2->length = 2;
  cf2->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 2);
  cf2->codedBases[0] = 0x92449244;
  cf2->codedBases[1] = 0x92449244;
  addChromToGenome(cf2, gf);

  strcpy(cf3->info, ">chr3");
  cf3->length = 3;
  cf3->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 3);
  cf3->codedBases[0] = 0x92449244;
  cf3->codedBases[1] = 0x92449244;
  cf3->codedBases[2] = 0x92449244;
  addChromToGenome(cf3, gf);

  strcpy(cf4->info, ">chr4");
  cf4->length = 4;
  cf4->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 4);
  cf4->codedBases[0] = 0x92449244;
  cf4->codedBases[1] = 0x92449244;
  cf4->codedBases[2] = 0x92449244;
  cf4->codedBases[3] = 0x92449244;
  addChromToGenome(cf4, gf);

  strcpy(cf5->info, ">chr5");
  cf5->length = 5;
  cf5->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 5);
  cf5->codedBases[0] = 0x92449244;
  cf5->codedBases[1] = 0x92449244;
  cf5->codedBases[2] = 0x92449244;
  cf5->codedBases[3] = 0x92449244;
  cf5->codedBases[4] = 0x92449244;
  addChromToGenome(cf5, gf);

  destroy_GenomeFa(gf);
}

static void _test_CodingBases() {
  char *bases[] = {"AAAACCCCGGGGTTTTAAAAC", "ACGTACGTACGTACGTACGTA", "ACGT"};
  uint64_t codedBases[] = {0x2494926DB9242494, 0x29C29C29C29C29C2,
                           0x29C0000000000000};
  ChromFa *cf = init_ChromFa();
  cf->codedBases = (uint64_t *)malloc(sizeof(uint64_t) * 3);
  cf->length = 46;
  strcpy(cf->info, ">codingBasesTestCf");

  // test the calculation of coded bases
  for (int i = 0; i < (sizeof(codedBases) / sizeof(uint64_t)); i++) {
    assert(codedBases[i] == codeBpBuf(bases[i]));
    cf->codedBases[i] = codedBases[i];
  }

  // test the extraction of coded bases from ChromFa
  assert_bases_equal(getBase(cf, 1), BASE_A);
  assert_bases_equal(getBase(cf, 10), BASE_G);
  assert_bases_equal(getBase(cf, 21), BASE_C);
  assert_bases_equal(getBase(cf, 22), BASE_A);
  assert_bases_equal(getBase(cf, 24), BASE_G);
  assert_bases_equal(getBase(cf, 42), BASE_A);
  assert_bases_equal(getBase(cf, 43), BASE_A);
  assert_bases_equal(getBase(cf, 44), BASE_C);
  assert_bases_equal(getBase(cf, 46), BASE_T);
  assert_bases_equal(getBase(cf, 999), BASE_INVALID);
}

static void _test_Loader() {
  FILE *faFile = fopen("data/testFa.fa", "r");
  if (faFile == NULL) {
    fprintf(stderr, "Error: failed to open file. \n");
    exit(EXIT_FAILURE);
  }
  GenomeFa *gf = init_GenomeFa();

  loadGenomeFaFromFile(gf, faFile);

  uint32_t testFaCfs_length[] = {151, 203, 152, 142};
  char *testFaCfs_info[] = {
      ">ch1  LN:151  AS:GRCh38", ">chr2  LN:203  AS:GRCh38",
      ">chr13  LN:152  AS:GRCh38", ">chr14  LN:142  AS:GRCh38"};

  ChromFa *tmpCf = gf->chroms;
  tmpCf = tmpCf->nextChrom;
  for (int i = 0; i < sizeof(testFaCfs_length) / sizeof(int); i++) {
    // printf("tmpCfLN: %d, testLN: %d\n", tmpCf->length, testFaCfs_length[i]);
    // printf("tmpCfinfo: \"%s\", testinfo: \"%s\"\n", tmpCf->info,
    // testFaCfs_info[i]);
    assert(tmpCf->length == testFaCfs_length[i]);
    assert(strcmp(tmpCf->info, testFaCfs_info[i]) == 0);
    tmpCf = tmpCf->nextChrom;
  }

  // TODO check the codedBases.
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
  tmpCf = gf->chroms->nextChrom;
  for (int i = 0; i < gf->chromNum; i++) {
    for (int j = 0; j < (tmpCf->length - 1) / BP_PER_UINT64; j++) {
      assert(tmpCf->codedBases[j] == codedChs[i][j]);
    }
    tmpCf = tmpCf->nextChrom;
  }

  destroy_GenomeFa(gf);
  fclose(faFile);
}

static void _test_Writer() {
  FILE *pFa = fopen("data/testFa.fa", "r");
  FILE *pO = fopen("data/output.fa", "w");  // It's an "O", not a "0"
  if (pFa == NULL || pO == NULL) {
    fprintf(stderr, "Error: failed to open file. \n");
    exit(EXIT_FAILURE);
  }

  GenomeFa *gf = init_GenomeFa();

  loadGenomeFaFromFile(gf, pFa);
  writeGenomeFaIntoFile(gf, pO);

  destroy_GenomeFa(gf);
  fclose(pFa);
  fclose(pO);
}

void _testSet_genomeFa() {
  _test_StructureLinks();
  _test_CodingBases();
  _test_Loader();
  _test_Writer();
}

// *********************** above are tests *************************

ChromFa *init_ChromFa() {
  ChromFa *cf = (ChromFa *)malloc(sizeof(ChromFa));
  if (cf == NULL) {
    fprintf(stderr, "Error: no enough memory for a new chrom.\n");
    exit(EXIT_FAILURE);
  }
  cf->length = 0;
  cf->codedBases = NULL;
  strcpy(cf->info, "initialized");
  cf->prevChrom = NULL;
  cf->nextChrom = NULL;
  return cf;
}

GenomeFa *init_GenomeFa() {
  GenomeFa *gf = (GenomeFa *)malloc(sizeof(GenomeFa));
  if (gf == NULL) {
    fprintf(stderr, "Error: no enough memory for a new genome.\n");
    exit(EXIT_FAILURE);
  }
  gf->chromNum = 0;
  gf->chroms = init_ChromFa();
  strcpy(gf->chroms->info, "header of genome fa");
  return gf;
}

ChromFa *getChromFromGenome_info(char *info, GenomeFa *gf) {
  if (gf == NULL) {
    fprintf(stderr, "Error: null pointer occurred for GenomeFa\n");
    exit(EXIT_FAILURE);
  }
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    if (strcmp(tmpCf->info, info) == 0) {
      return tmpCf;
    }
    tmpCf = tmpCf->nextChrom;
  }
  return NULL;
}

ChromFa *getChromFromGenome_Idx(int idx, GenomeFa *gf) {
  // TODO not tested
  if (gf == NULL) {
    fprintf(stderr, "Error: null pointer occurred for GenomeFa\n");
    exit(EXIT_FAILURE);
  }
  ChromFa *tmpCf = gf->chroms;
  int tmpIdx = 0;
  while (tmpCf != NULL) {
    if (tmpIdx == idx) {
      return tmpCf;
    }
    tmpIdx++;
    tmpCf = tmpCf->nextChrom;
  }

  return NULL;
}

Base getBase(ChromFa *cf, uint32_t pos) {
  static uint64_t uint64Length = sizeof(uint64_t);
  uint32_t arrayLength = cf->length / uint64Length;
  uint32_t arrayIdx = pos / (BP_PER_UINT64 + 1);

  if (arrayIdx > arrayLength || pos <= 0) {
    return BASE_INVALID;
  } else {
    // The following commented lines were used for debugging.
    // uint64_t arrayElement = cf->codedBases[arrayIdx];
    // // 1 % 21 = 1, 21 % 21 = 0. This is not good. So (21 - 1) % 21 + 1 = 21.
    // uint8_t baseIdx = (pos - 1) % BP_PER_UINT64 + 1;
    // uint8_t offset = 1 + BASE_CODE_LENGTH * (BP_PER_UINT64 - baseIdx);

    // Base retVal = (Base)((arrayElement >> offset) & BASE_MASK_RIGHT);
    return (Base)(
        (cf->codedBases[arrayIdx] >>
             (1 + BASE_CODE_LENGTH *
                      (BP_PER_UINT64 - ((pos - 1) % BP_PER_UINT64 + 1))) &
         BASE_MASK_RIGHT));
  }
}

void addChromToGenome(ChromFa *cf, GenomeFa *gf) {
  if (cf == NULL || gf == NULL) {
    fprintf(stderr, "Error: null pointer occurred for ChromFa or GenomeFa\n");
    exit(EXIT_FAILURE);
  }
  if (gf->chroms == NULL) {
    gf->chroms = cf;
    gf->chromNum++;
    return;
  } else {
    ChromFa *tmpCf = gf->chroms;
    while (tmpCf->nextChrom != NULL) {
      tmpCf = tmpCf->nextChrom;
    }
    tmpCf->nextChrom = cf;
    cf->prevChrom = tmpCf;
  }
}

void destroy_ChromFa(ChromFa *cf) {
  if (cf == NULL) {
    fprintf(stderr,
            "Error: null pointer occurred when destroying a ChromFa object\n");
    exit(EXIT_FAILURE);
  }
  ChromFa *prevCf = cf->prevChrom;
  ChromFa *nextCf = cf->nextChrom;
  if (prevCf != NULL) {
    prevCf->nextChrom = nextCf;
  }
  if (nextCf != NULL) {
    nextCf->prevChrom = prevCf;
  }
  free(cf->codedBases);
  free(cf);
}

void destroy_GenomeFa(GenomeFa *gf) {
  if (gf == NULL) {
    fprintf(stderr,
            "Error: null pointer occurred when destroying a GenomeFa object\n");
    exit(EXIT_FAILURE);
  }
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    ChromFa *nextCf = tmpCf->nextChrom;
    destroy_ChromFa(tmpCf);
    tmpCf = nextCf;
  }
  free(gf);
}

void printGenomeFa(GenomeFa *gf) {
  printf("number of chroms: %d\n", gf->chromNum);
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    printChromFa(tmpCf);
    tmpCf = tmpCf->nextChrom;
  }
}

void printGenomeFa_brief(GenomeFa *gf) {
  ChromFa *tmpCf = gf->chroms;
  while (tmpCf != NULL) {
    printf("%s\n", tmpCf->info);
    tmpCf = tmpCf->nextChrom;
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

uint64_t codeBpBuf(char *bpBuf) {
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

int newInfoForGenomeFa(GenomeFa *gf, char *infoBuf) {
  ChromFa *cf = NULL;
  if (getChromFromGenome_info(infoBuf, gf) != NULL) {
    fprintf(stderr, "Error: duplicate chromosomes. \n");
    exit(EXIT_FAILURE);
  } else {
    cf = init_ChromFa();
  }

  // Note that the following lines are coded according to the format of standard
  // *.fa files, or to say, the GRCh38.fa file. The main purpose is to extract
  // the "LN" field from the info line.
  int infoLen = strlen(infoBuf);
  if (infoLen > 0 && infoBuf[0] != '>') {
    return 0;
  }

  int chromLN = 0;
  static char numberBuf[40];
  int takeInNumberFlag = 0;
  int numberBufIdx = 0;

  for (int i = 0; i < infoLen; i++) {
    if ((i + 2) < infoLen && infoBuf[i] == 'L' && infoBuf[i + 1] == 'N' &&
        infoBuf[i + 2] == ':') {
      // jump to the number part if it is the "LN:..."
      i = i + 3;
      takeInNumberFlag = 1;
    }
    if (takeInNumberFlag == 1) {
      if (infoBuf[i] >= '0' && infoBuf[i] <= '9') {
        numberBuf[numberBufIdx] = infoBuf[i];
        numberBufIdx++;
      } else {
        numberBuf[numberBufIdx] = '\0';
        chromLN = atoi(numberBuf);
        cf->length = chromLN;
        strcpy(cf->info, infoBuf);
        cf->codedBases = (uint64_t *)malloc(
            sizeof(uint64_t) * ((chromLN - 1) / BP_PER_UINT64 + 1));
        addChromToGenome(cf, gf);
        return 1;
      }
    }
  }

  return 0;
}

void loadGenomeFaFromFile(GenomeFa *gf, FILE *fp) {
  if (fp == NULL) {
    fprintf(stderr, "Error: empty file pointer. \n");
    exit(EXIT_FAILURE);
  }
  if (gf == NULL) {
    fprintf(stderr, "Error: null pointer for GenomeFa\n");
    exit(EXIT_FAILURE);
  }

  char bpBuf[BP_PER_UINT64 + 1];    // buffer for bases in *.fa file
  char infoBuf[MAX_RECORD_LENGTH];  // buffer for info line in *.fa file
  int bpBufIdx = 0;
  int infoBufIdx = 0;
  int bpPosWithinChrom =
      1;  // 1-base position of bases within chrom. The index of coded bases
          // element index in the ChromFa object needs to be calculated
  int isBpLine = 0;  // if the current line contains the bases
  int ifIgnore = 0;  // selectively ignore chroms whose info lines do not match
                     // the specifications
  ChromFa *currentCf = NULL;  // current chrom to be loaded from file. This can
                              // accelerate the program by cutting off repeative
                              // steps of getting the ChromFa object
  char tmpCh;
  clock_t start = clock(), end = 0;
  while ((tmpCh = fgetc(fp)) != EOF) {
    // use some little tricks to set the value of isBpLine and initialize the
    // GenomeFa object for following base-coding process
    if (tmpCh == '>') {     // reach the start of the info line
      if (bpBufIdx != 0) {  // flush the uncoded bases in the bpBuf. This is
                            // executed when the bpBuf is not filled up but the
                            // chromosome has already reached its end.
        bpBuf[bpBufIdx] = '\0';
        currentCf->codedBases[(bpPosWithinChrom - 1) / BP_PER_UINT64] =
            codeBpBuf(bpBuf);
        bpBufIdx = 0;
      }
      isBpLine = 0;
      infoBufIdx = 0;
      infoBuf[infoBufIdx++] = tmpCh;
      ifIgnore = 0;
      continue;
    }
    if (ifIgnore == 1) {
      continue;
    }
    if (isBpLine == 0) {
      if (tmpCh == '\n') {           // reach the end of the info line
        infoBuf[infoBufIdx] = '\0';  // pad the end of a string
        // pass the info string to the helper function for process
        if (newInfoForGenomeFa(gf, infoBuf) == 0) {
          // if the format does not match the specifications, ignore it
          fprintf(stderr, "chrom ignored: \"%s\"\n", infoBuf);
          ifIgnore = 1;
          continue;
        }
        gf->chromNum++;
        if ((currentCf = getChromFromGenome_info(infoBuf, gf)) == NULL) {
          fprintf(stderr,
                  "Error: failed to locate the chrom with info \'%s\'\n",
                  infoBuf);
          exit(EXIT_FAILURE);
        }
        isBpLine = 1;
        bpPosWithinChrom =
            1;  // forgetting about this will result in an error when malloc()
                // memory for a new ChromFa (I don't know why)
        continue;
      } else {  // extract info line from *.fa file and put into infoBuf
        infoBuf[infoBufIdx++] = tmpCh;
        continue;
      }
    }  // end of (isBpLine == 0)

    // code the bases into binary integers and store them into the
    // designated data structure properly.
    if (isBpLine == 1) {
      if (tmpCh == '\n') {
        continue;
      } else {
        bpBuf[bpBufIdx++] = tmpCh;
        if (bpBufIdx == BP_PER_UINT64) {  // when the bpBuf is full
          bpBuf[bpBufIdx] = '\0';
          currentCf->codedBases[(bpPosWithinChrom - 1) / BP_PER_UINT64] =
              codeBpBuf(bpBuf);
          bpBufIdx = 0;
        }
        bpPosWithinChrom++;
      }
    }  // end of (isBpLine == 1)
  }
  // flush the uncoded bases in the bpBuf. This is executed when the bpBuf is
  // not filled up but the chromosome has already reached its end.
  if (bpBufIdx != 0) {
    bpBuf[bpBufIdx] = '\0';
    currentCf->codedBases[(bpPosWithinChrom - 1) / BP_PER_UINT64] =
        codeBpBuf(bpBuf);
    bpBufIdx = 0;
  }
  end = clock();
  assert(printf("... fa data loading finished. total time(s): %f\n",
         (double)(end - start) / CLOCKS_PER_SEC)>=0);
}

void writeGenomeFaIntoFile(GenomeFa *gf, FILE *fp) {
  // TODO
  ChromFa *tmpCf = gf->chroms->nextChrom;

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
    tmpCf = tmpCf->nextChrom;
  }
}
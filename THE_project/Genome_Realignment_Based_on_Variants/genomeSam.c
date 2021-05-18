#include "genomeSam.h"

char *rsDataSeq(RecSam *rs) {
  uint32_t seqLength = rs->rec->core.l_qseq;
  char *seq = (char *)calloc(seqLength + 1, sizeof(char));
  for (int i = 0; i < seqLength; i++) {
    char bp;
    switch (bam_seqi(bam_get_seq(rsData(rs)), i)) {
      case 1: {
        bp = 'A';
        break;
      }
      case 2: {
        bp = 'C';
        break;
      }
      case 4: {
        bp = 'G';
        break;
      }
      case 8: {
        bp = 'T';
        break;
      }
      case 15: {
        bp = 'N';
        break;
      }
    }
    seq[i] = bp;
  }
  seq[seqLength] = '\0';
  return seq;
}

static RecSam *init_RecSam() {
  RecSam *rs = (RecSam *)malloc(sizeof(RecSam));
  rs->rec = NULL;
  rs->next = NULL;
  return rs;
}

static void destroy_RecSam(RecSam *rs) {
  if (rs == NULL) return;
  bam_destroy1(rs->rec);
  free(rs);
}

static ChromSam *init_ChromSam() {
  ChromSam *cs = (ChromSam *)malloc(sizeof(ChromSam));
  cs->name = "";
  cs->recCnt = 0;
  cs->rss = init_RecSam();
  cs->rss->rec = bam_init1();
  cs->rss->rec->core.pos = -2;  // smaller than any record's pos in a vcf file
  cs->next = NULL;
  return cs;
}

static void destroy_ChromSam(ChromSam *cs) {
  if (cs == NULL) return;
  RecSam *rs = cs->rss;
  RecSam *tmpRs = NULL;
  uint32_t hit = 0;
  while (rs != NULL) {
    hit++;
    tmpRs = rs->next;
    destroy_RecSam(rs);
    cs->recCnt--;
    rs = tmpRs;
  }
  free(cs->name);
  free(cs);
}

GenomeSam *init_GenomeSam() {
  GenomeSam *gs = (GenomeSam *)malloc(sizeof(GenomeSam));
  if (gs == NULL) {
    fprintf(stderr, "Error: memory not enough for new GenomeSam object. \n");
    exit(EXIT_FAILURE);
  }
  gs->chromCnt = 0;
  gs->hdr = NULL;
  gs->css = NULL;
  return gs;
}

void destroy_GenomeSam(GenomeSam *gs) {
  if (gs == NULL) return;
  sam_hdr_destroy(gs->hdr);

  ChromSam *cs = gs->css;
  while (cs != NULL) {
    ChromSam *tmpCs = cs->next;
    destroy_ChromSam(cs);
    gs->chromCnt--;
    cs = tmpCs;
  }
  free(gs);
}

GenomeSamIterator *init_GenomeSamIterator(GenomeSam *gs) {
  GenomeSamIterator *gsIt =
      (GenomeSamIterator *)malloc(sizeof(GenomeSamIterator));
  gsIt->gs = gs;
  gsIt->tmpCs = NULL;
  gsIt->tmpRs = NULL;
  return gsIt;
}

void destroy_GenomeSamIterator(GenomeSamIterator *gsIt) { free(gsIt); }

ChromSam *gsItNextChrom(GenomeSamIterator *gsIt) {
  if (gsIt->gs == NULL) {
    fprintf(stderr, "Warning: gsIterator not initailized properly. \n");
    return NULL;
  }
  if (gsIt->tmpCs != NULL) {
    gsIt->tmpCs = gsIt->tmpCs->next;
  } else {
    gsIt->tmpCs = gsIt->gs->css;
  }
  return gsIt->tmpCs;
}

RecSam *gsItNextRec(GenomeSamIterator *gsIt) {
  if (gsIt->gs == NULL) {
    fprintf(stderr, "Warning: gsIterator not initailized properly. \n");
    return NULL;
  }
  if (gsIt->tmpRs != NULL) {
    gsIt->tmpRs = gsIt->tmpRs->next;
  } else {
    if (gsIt->tmpCs != NULL) {
      gsIt->tmpRs = gsIt->tmpCs->rss;
      // Remember that the rss in cs is a linked-list with an empty header
      gsIt->tmpRs = gsIt->tmpRs->next;
    } else {
      gsIt->tmpRs = NULL;
    }
  }
  return gsIt->tmpRs;
}

void addChromToGenomeSam(ChromSam *cs, GenomeSam *gs) {
  ChromSam *tmpCs = gs->css;
  if (tmpCs == NULL) {  // if there is no chrom in GenomeSam
    gs->css = cs;
    gs->chromCnt++;
    return;
  }

  while (tmpCs->next != NULL) {
    if (strcmp(cs->name, tmpCs->name) ==
        0) {  // if there already exists the same cs
      fprintf(stderr, "Warning: trying to add duplicated ChromSam. \n");
      return;
    }
    tmpCs = tmpCs->next;
  }
  // if same cs not found, add it to the end of the linked-list of ChromSam
  tmpCs->next = cs;
  gs->chromCnt++;
}

void addRecToChromSam(RecSam *rs, ChromSam *cs) {
  RecSam *tmpRec = cs->rss->next;
  RecSam *lastRec = cs->rss;

  if (tmpRec == NULL) {
    lastRec->next = rs;
    cs->recCnt++;
    return;
  }
  while (tmpRec != NULL) {
    if (getRecSam_pos(tmpRec) <= getRecSam_pos(rs)) {
      lastRec = tmpRec;
      tmpRec = tmpRec->next;
    } else {  // when pos(rs) < pos(tmpRec), insert rs between
              // "...,lastRec,tmpRec,..."
      lastRec->next = rs;
      rs->next = tmpRec;
      cs->recCnt++;
      return;
    }
  }
  // if always pos(tmpRec) <= pos(rs), add rs to the end
  lastRec->next = rs;
  cs->recCnt++;
}

ChromSam *getChromFromGenomeSam(char *chromName, GenomeSam *gs) {
  ChromSam *tmpCs = gs->css;
  while (tmpCs != NULL) {
    if (strcmp(tmpCs->name, chromName) == 0) {
      return tmpCs;
    } else {
      tmpCs = tmpCs->next;
    }
  }
  // if never found ChromSam with the same name as chromName
  return NULL;
}

RecSam *getRecFromChromSam(uint32_t idx, ChromSam *cs) {
  // TODO not tested
  RecSam *tmpRs = cs->rss->next;  // pass the header rv
  if (idx >= cs->recCnt) {
    fprintf(stderr,
            "Error: index out of bounder when getting record from a chrom. \n");
    exit(EXIT_FAILURE);
  }
  uint32_t tmpIdx = 0;
  while (tmpRs != NULL) {
    if (tmpIdx == idx) {
      return tmpRs;
    } else {
      tmpIdx++;
    }
  }
}

char *getRecSam_chNam(RecSam *rs, GenomeSam *gs) {
  const char *chNam = sam_hdr_tid2name(gs->hdr, rs->rec->core.tid);
  if (chNam == NULL) {
    /*
     * This part of code was originally only "return 'unknown'". And that
     * resulted in a bug. Because you cannot free a string that is pre-allocated
     * by the compiler instead allocated dynamically.
     */
    static const int len = strlen("unknown)");
    char *unknownNam = (char *)malloc(sizeof(char) * (len + 1));
    strcpy(unknownNam, "(unknown)");
    return unknownNam;
  } else {
    return strdup(chNam);
  }
}

void loadGenomeSamFromFile(GenomeSam *gs, char *filePath) {
  htsFile *fp = hts_open(filePath, "r");
  sam_hdr_t *hdr = sam_hdr_read(fp);
  bam1_t *rec = bam_init1();

  if (hdr == NULL) {
    fprintf(stderr, "Error: failed creating bcf header struct.\n");
    exit(EXIT_FAILURE);
  } else {
    gs->hdr = sam_hdr_dup(hdr);
  }
  if (rec == NULL) {
    fprintf(stderr, "Error: memory not enough for new bcf1_t object.\n");
    exit(EXIT_FAILURE);
  }

  ChromSam *lastUsedChrom = NULL;
  uint32_t loadedCnt = 0;
  while (sam_read1(fp, hdr, rec) >= 0) {
    // printSamRecord_brief(hdr, rec);
    loadedCnt++;
    RecSam *newRs = init_RecSam();
    newRs->rec = bam_dup1(rec);

    //  add the record into GenomeSam object
    char *rsChNam = getRecSam_chNam(newRs, gs);
    if (lastUsedChrom != NULL &&
        strcmp(rsChNam, lastUsedChrom->name) ==
            0) {  // if the new record points to the same chrom as the last one
      addRecToChromSam(newRs, lastUsedChrom);
    } else {
      // if the new record points to another chrom
      lastUsedChrom = getChromFromGenomeSam(rsChNam, gs);
      if (lastUsedChrom == NULL) {
        // if there is no such chrom as the new record points to
        ChromSam *newCs = init_ChromSam();
        newCs->name = rsChNam;
        addChromToGenomeSam(newCs, gs);
        lastUsedChrom = newCs;
        // printf("... new chrom assigned, name: %s\n", rsChNam);
      }
      addRecToChromSam(newRs, lastUsedChrom);
    }
  }

  bam_destroy1(rec);
  sam_hdr_destroy(hdr);
  hts_close(fp);
}

void writeGenomeSamIntoFile(GenomeSam *gs, char *filePath) {
  htsFile *fp = hts_open(filePath, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filePath);
    exit(EXIT_FAILURE);
  }
  GenomeSamIterator *gsIt = init_GenomeSamIterator(gs);
  ChromSam *tmpCs = gsItNextChrom(gsIt);
  RecSam *tmpRs = gsItNextRec(gsIt);

  // write file headers and macros
  if(sam_hdr_write(fp, gsDataHdr(gs)) < 0){
    fprintf(stderr, "Error: failed to write sam file header. \n");
    exit(EXIT_FAILURE);
  }

  // write sam records
  while (tmpRs != NULL) {
    // printSamRecord_brief(gs, tmpRs->rec);
    bam1_t *rec = rsData(tmpRs);
    if (sam_write1(fp, gsDataHdr(gs), rec) < 0) {
      fprintf(stderr, "Error: failed to write sam record.\n");
      exit(EXIT_FAILURE);
    }

    tmpRs = gsItNextRec(gsIt);
    if (tmpRs == NULL) {
      tmpCs = gsItNextChrom(gsIt);
      tmpRs = gsItNextRec(gsIt);
    }
  }

  destroy_GenomeSamIterator(gsIt);
  hts_close(fp);
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

static int _test_LoadingAndIterator() {
  GenomeSam *gs = init_GenomeSam();

  loadGenomeSamFromFile(gs, "data/example.sam");

  // printGenomeSam(gs);

  GenomeSamIterator *gsIt = init_GenomeSamIterator(gs);
  ChromSam *tmpCs = gsItNextChrom(gsIt);
  RecSam *tmpRs = gsItNextRec(gsIt);

  while (tmpRs != NULL) {
    // printSamRecord_brief(gs, tmpRs->rec);

    tmpRs = gsItNextRec(gsIt);
    if (tmpRs == NULL) {
      tmpCs = gsItNextChrom(gsIt);
      tmpRs = gsItNextRec(gsIt);
    }
  }

  destroy_GenomeSamIterator(gsIt);
  destroy_GenomeSam(gs);
  return 1;
}

static int _test_WritingAndIterator() {
  char *outputFilePath = "data/exampleOut.sam";
  GenomeSam *gs = init_GenomeSam();

  loadGenomeSamFromFile(gs, "data/example.sam");

  writeGenomeSamIntoFile(gs, outputFilePath);

  destroy_GenomeSam(gs);
  return 1;
}

void _testSet_genomeSam() {
  assert(_test_LoadingAndIterator());
  assert(_test_WritingAndIterator());
}

void printSamHeader(bam_hdr_t *header) {
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

void printSamRecord_brief(GenomeSam *gs, bam1_t *record) {
  bam_hdr_t *hdr = gs->hdr;
  // qname flag rname pos mapq
  printf("%s\t", bam_get_qname(record));
  // flag
  printf("0x%" PRIx16 "\t", record->core.flag);
  // rname
  printf("%s\t", sam_hdr_tid2name(hdr, record->core.tid));
  // pos: 1-based
  printf("%" PRId64 "\t", record->core.pos + 1);
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
  // some differences on the definition of this field. Ignore it
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
  printf("\n");
}

void printGenomeSam(GenomeSam *gs) {
  if (gs == NULL) return;
  printSamHeader(gs->hdr);
  ChromSam *tmpCs = gs->css;
  while (tmpCs != NULL) {
    printChromSam(gs, tmpCs);
    tmpCs = tmpCs->next;
  }
}

void printChromSam(GenomeSam *gs, ChromSam *cs) {
  if (cs == NULL || gs == NULL) return;
  printf("chrom: %s, recCnt: %" PRIu32 "\n", cs->name, cs->recCnt);
  RecSam *rs = cs->rss->next;
  while (rs != NULL) {
    printSamRecord_brief(gs, rsData(rs));
    rs = rs->next;
  }
}
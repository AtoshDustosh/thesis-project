#include "genomeVcf.h"

void _testSet_genomeVcf() {
  GenomeVcf *gv = init_GenomeVcf();

  loadGenomeVcfFromFile(gv, "data/example.vcf");

  // printGenomeVcf(gv);

  destroy_GenomeVcf(gv);
}

void printVcfHeader(bcf_hdr_t *hdr) {
  printf("vcf version: %s\n", bcf_hdr_get_version(hdr));
  printf("\n");
}

void printVcfRecord(bcf1_t *rec) {
  /*
   * The ids of chroms correspond to the sequence they appear in the
   * vcf/bcf file. For example, if contigs are written as follows,
   * ##contig=<ID=5>
   * ##contig=<ID=1>
   * ##contig=<ID=3>
   * id array of the 3 chroms/contigs should be 0, 1, 2, instead of
   * 3, 1, 2.
   */
  printf("chrom: %d\t", rec->rid);
  /*
   * The output on console looks like 0-based index. But that's
   * for the convenience of array in C.
   */
  printf("pos: %ld\t", rec->pos);
  printf("qual: %f\n", rec->qual);
  /*
   * Including the reference allele, which is "record->d.allele[0]".
   */
  printf("number of alleles: %d\n", rec->n_allele);
  printf("alleles: ");
  for (int i = 0; i < rec->n_allele; i++) {
    printf("%s ", rec->d.allele[i]);
  }
  printf("\n");

  printf("variant types:");
  for (int i = 1; i < rec->n_allele; i++) {
    printf("%d ", bcf_get_variant_type(rec, i));
  }
  printf("\n");
  // other fields are temporarily not needed
  printf("\n");
}

void printVcfRecord_brief(bcf_hdr_t *hdr, bcf1_t *rec) {
  // chrom
  printf("%s\t", bcf_seqname_safe(hdr, rec));
  // pos
  printf("%" PRId64 "\t", rec->pos);
  // id
  printf("%s\t", rec->d.id);
  // ref alt
  printf("%s\t", rec->d.allele[0]);
  if (rec->n_allele == 1) {
    printf(".");
  } else {
    for (int i = 1; i < rec->n_allele; i++) {
      if (i >= 2) {
        printf(",%s", rec->d.allele[i]);
      } else {
        printf("%s", rec->d.allele[i]);
      }
    }
  }
  printf("\t");
  // qual
  printf("%f\t", rec->qual);
  // filter, info, format and other fields are ignored
  printf("\n");
}

void printGenomeVcf(GenomeVcf *gv) {
  if (gv == NULL) return;
  printVcfHeader(gv->hdr);
  ChromVcf *tmpCv = gv->cv;
  while (tmpCv != NULL) {
    printChromVcf(gv, tmpCv);
    tmpCv = tmpCv->next;
  }
}

void printChromVcf(GenomeVcf *gv, ChromVcf *cv) {
  if (cv == NULL || gv == NULL) return;
  printf("chrom: %s, recCnt: %" PRIu32 "\n", cv->name, cv->recCnt);
  RecVcf *rv = cv->rv->next;
  while (rv != NULL) {
    printVcfRecord_brief(gv->hdr, rv->rec);
    rv = rv->next;
  }
}

static RecVcf *init_RecVcf() {
  RecVcf *rv = (RecVcf *)malloc(sizeof(RecVcf));
  rv->rec = NULL;
  rv->next = NULL;
  return rv;
}

static void destroy_RecVcf(RecVcf *rv) {
  if (rv == NULL) return;
  bcf_destroy1(rv->rec);
  free(rv);
}

static ChromVcf *init_ChromVcf() {
  ChromVcf *cv = (ChromVcf *)malloc(sizeof(ChromVcf));
  cv->name = "";
  cv->recCnt = 0;
  cv->rv = init_RecVcf();
  cv->rv->rec = bcf_init1();
  cv->rv->rec->pos = -1;  // smaller than any record's pos in a vcf file
  cv->next = NULL;
  return cv;
}

static void destroy_ChromVcf(ChromVcf *cv) {
  if (cv == NULL) return;
  RecVcf *rv = cv->rv;
  while (rv != NULL) {
    RecVcf *tmpRv = rv->next;
    destroy_RecVcf(rv);
    cv->recCnt--;
    rv = tmpRv;
  }
  free(cv->name);
  free(cv);
}

GenomeVcf *init_GenomeVcf() {
  GenomeVcf *gv = (GenomeVcf *)malloc(sizeof(GenomeVcf));
  if (gv == NULL) {
    fprintf(stderr, "Error: memory not enough for new GenomeVcf object. \n");
    exit(EXIT_FAILURE);
  }
  gv->chromCnt = 0;
  gv->hdr = NULL;
  gv->cv = NULL;
  return gv;
}

void destroy_GenomeVcf(GenomeVcf *gv) {
  if (gv == NULL) return;
  bcf_hdr_destroy(gv->hdr);

  ChromVcf *cv = gv->cv;
  while (cv != NULL) {
    ChromVcf *tmpCv = cv->next;
    destroy_ChromVcf(cv);
    gv->chromCnt--;
    cv = tmpCv;
  }
  free(gv);
}

void addChromToGenomeVcf(ChromVcf *cv, GenomeVcf *gv) {
  ChromVcf *tmpCv = gv->cv;
  if (tmpCv == NULL) {  // if there is no chrom in GenomeVcf
    gv->cv = cv;
    gv->chromCnt++;
    return;
  }

  while (tmpCv->next != NULL) {
    if (strcmp(cv->name, tmpCv->name) ==
        0) {  // if there already exists the same cv
      fprintf(stderr, "Warning: trying to add duplicated ChromVcf. \n");
      return;
    }
    tmpCv = tmpCv->next;
  }
  // if same cv not found, add it to the end of the linked-list of ChromVcf
  tmpCv->next = cv;
  gv->chromCnt++;
}

void addRecToChromVcf(RecVcf *rv, ChromVcf *cv) {
  RecVcf *tmpRec = cv->rv->next;
  RecVcf *lastRec = cv->rv;

  if (tmpRec == NULL) {
    lastRec->next = rv;
    cv->recCnt++;
    return;
  }
  while (tmpRec != NULL) {
    if (getRecVcf_pos(tmpRec) <= getRecVcf_pos(rv)) {
      lastRec = tmpRec;
      tmpRec = tmpRec->next;
    } else {  // when pos(rv) < pos(tmpRec), insert rv between
              // "...,lastRec,tmpRec,..."
      lastRec->next = rv;
      rv->next = tmpRec;
      cv->recCnt++;
      return;
    }
  }
  // if always pos(tmpRec) <= pos(rv), add rv to the end
  lastRec->next = rv;
  cv->recCnt++;
}

ChromVcf *getChromFromGenomeVcf(char *chromName, GenomeVcf *gv) {
  ChromVcf *tmpCv = gv->cv;
  while (tmpCv != NULL) {
    if (strcmp(tmpCv->name, chromName) == 0) {
      return tmpCv;
    } else {
      tmpCv = tmpCv->next;
    }
  }
  // if never found ChromVcf with the same name as chromName
  return NULL;
}

RecVcf *getRecFromChromVcf(uint32_t idx, ChromVcf *cv) {
  // TODO not tested
  RecVcf *tmpRv = cv->rv->next;  // pass the header rv
  if (idx >= cv->recCnt) {
    fprintf(stderr,
            "Error: index out of bounder when getting record from a chrom. \n");
    exit(EXIT_FAILURE);
  }
  uint32_t tmpIdx = 0;
  while (tmpRv != NULL) {
    if (tmpIdx == idx) {
      return tmpRv;
    } else {
      tmpIdx++;
    }
  }
}

char *getRecVcf_chNam(RecVcf *rv, GenomeVcf *gv) {
  return strdup(bcf_seqname_safe(gv->hdr, rv->rec));
}

void loadGenomeVcfFromFile(GenomeVcf *gv, char *filePath) {
  htsFile *fp = hts_open(filePath, "r");
  bcf_hdr_t *hdr = bcf_hdr_read(fp);
  bcf1_t *rec = bcf_init1();

  if (hdr == NULL) {
    fprintf(stderr, "Error: failed creating bcf header struct.\n");
    exit(EXIT_FAILURE);
  } else {
    gv->hdr = bcf_hdr_dup(hdr);
  }
  if (rec == NULL) {
    fprintf(stderr, "Error: memory not enough for new bcf1_t object.\n");
    exit(EXIT_FAILURE);
  }

  clock_t start = clock(), end = 0;
  ChromVcf *lastUsedChrom = NULL;
  uint32_t loadedCnt = 0;
  printf("Loading %s ... \n", filePath);
  while (bcf_read1(fp, hdr, rec) >= 0) {
    // The "bcf_unpack" method must be called for every new bcf1_t object
    bcf_unpack(rec, BCF_UN_ALL);
    bcf1_t *tmpRec = bcf_dup(rec);
    bcf_unpack(tmpRec, BCF_UN_ALL);
    // printVcfRecord_brief(hdr, tmpRec);
    loadedCnt++;

    RecVcf *newRv = init_RecVcf();
    newRv->rec = tmpRec;
    // TODO add the record into GenomeVcf object
    char *rvChNam = getRecVcf_chNam(newRv, gv);
    if (lastUsedChrom != NULL &&
        strcmp(rvChNam, lastUsedChrom->name) ==
            0) {  // if the new record points to the same chrom as the last one
      addRecToChromVcf(newRv, lastUsedChrom);
    } else {
      // if the new record points to another chrom
      lastUsedChrom = getChromFromGenomeVcf(rvChNam, gv);
      if (lastUsedChrom == NULL) {
        // if there is no such chrom as the new record points to
        ChromVcf *newCv = init_ChromVcf();
        newCv->name = rvChNam;
        addChromToGenomeVcf(newCv, gv);
        lastUsedChrom = newCv;
        // printf("... new chrom assigned, name: %s\n", rvChNam);
      }
      addRecToChromVcf(newRv, lastUsedChrom);
    }
  }
  end = clock();
  printf("... vcf data loading finished. Processed %" PRIu32
         " records. Total time(s): %f\n",
         loadedCnt, (double)(end - start) / CLOCKS_PER_SEC);

  bcf_destroy1(rec);
  bcf_hdr_destroy(hdr);
  hts_close(fp);
}

void writeGenomeVcfIntoFile(GenomeVcf *gv, char *filePath) {
  // TODO
  fprintf(stderr, "Warning: method writeGenomeVcfIntoFile not implemented. \n");
}
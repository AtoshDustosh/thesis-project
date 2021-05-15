#include "integrateVcfToSam.h"

void check_files(Options *opts) {
  if (getSamFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete for variants integration - lak input "
            "files *.sam/bam\n");
    exit(EXIT_FAILURE);
  }
  if (getVcfFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete for variants integration - lak input "
            "files *.vcf/bcf\n");
    exit(EXIT_FAILURE);
  }
  if (getFaFile(opts) == NULL) {
    fprintf(stderr,
            "Error: arguments incomplete for variants integration - lak input "
            "files *.fa/fna/fasta\n");
    exit(EXIT_FAILURE);
  }
  if (getOutputFile(opts) == NULL) {
    fprintf(stderr,
            "Warning: output file not set. Set as default output file: "
            "%s\n",
            default_outputFile);
    setOutputFile(opts, default_outputFile);
  }
}

typedef struct _define_ThreadArgs {
  int64_t id;  // identifier for the thread
  Options *opts;
  GenomeFa *gf;
  GenomeSam *gs;
  GenomeVcf_bplus *gv;
  int64_t id_sam_start;  // index for sam record when the thread starts
  int64_t id_sam_end;    // index for sam record when the thread ends
} ThreadArgs;

void *integration_threads(void *args) {
  ThreadArgs *args_thread = (ThreadArgs *)args;

  // printf("thread (%" PRId64 ") process sam rec from %" PRId64 " to %" PRId64
  //        "\n",
  //        args_thread->id, args_thread->id_sam_start,
  //        args_thread->id_sam_end);
  // Execute integrations process

  char path_outputFile[strlen(getOutputFile(args_thread->opts)) + 64];
  sprintf(path_outputFile, "%s.thread%" PRId64 "",
          getOutputFile(args_thread->opts), args_thread->id);

  printf("output file name for thread (%" PRId64 "): %s\n", args_thread->id,
         path_outputFile);

  samFile *file_output = NULL;
  file_output = sam_open(getOutputFile(args_thread->opts), "w");
  if (file_output == NULL) {
    fprintf(stderr, "Error: cannot open file %s with mode \"w\"\n",
            path_outputFile);
    exit(EXIT_FAILURE);
  }

  // Write header into file
  sam_hdr_t *header = gsDataHdr(args_thread->gs);
  if (sam_hdr_write(file_output, header) < 0) {
    fprintf(stderr, "Error: failed writing header for %s\n", path_outputFile);
    exit(EXIT_FAILURE);
  }

  // Extract arguments from thread input
  GenomeFa *gf = args_thread->gf;
  GenomeSam *gs = args_thread->gs;
  GenomeVcf_bplus *gv = args_thread->gv;
  
  const int sv_min_len = getSVminLen(args_thread->opts);

  // Iterate all sam records and locate their corresponding variants.
  GenomeSamIterator *gsIt = init_GenomeSamIterator(gs);
  ChromSam *cs_tmp = gsItNextChrom(gsIt);
  RecSam *rs_tmp = gsItNextRec(gsIt);

  int64_t id_rec = 0;
  int64_t id_rec_start = args_thread->id_sam_start;  // included
  int64_t id_rec_end = args_thread->id_sam_end;      // included
  while (rs_tmp != NULL) {
    // -------- only handle records within [id_start, id_end] --------
    if (id_rec < id_rec_start) {
      rs_tmp = gsItNextRec(gsIt);
      if (rs_tmp == NULL) {
        cs_tmp = gsItNextChrom(gsIt);
        rs_tmp = gsItNextRec(gsIt);
      }
      id_rec++;
      continue;
    }
    if (id_rec > id_rec_end) {
      break;
    }
    // ---------- get information of temporary sam record ------------
    const char *rname_read = rsDataRname(gs, rs_tmp);
    const char *qname_read = rsDataQname(rs_tmp);
    int64_t pos_start_read = rsDataPos(rs_tmp);
    uint32_t length_read = rsDataSeqLength(rs_tmp);
    char *seq_read = rsDataSeq(rs_tmp);
    
    // TODO iterator for variants on the vcf bplus structure
    // TODO add field for RecVcf_bplus (bpnode*) to find its bpnode
    static int test_integer;

    printf("test_integer: %d\n", test_integer++);

    // ----------------------- free memories -------------------------
    free(seq_read);

    // --------------------- keep on iterating -----------------------
    rs_tmp = gsItNextRec(gsIt);
    if (rs_tmp == NULL) {
      cs_tmp = gsItNextChrom(gsIt);
      rs_tmp = gsItNextRec(gsIt);
    }
    id_rec++;
  }

  destroy_GenomeSamIterator(gsIt);

  sam_close(file_output);

  return (void *)(args_thread->id);
}

void integration(Options *opts) {
  check_files(opts);

  // Init structures (data storage and access)
  GenomeFa *gf = init_GenomeFa();
  loadGenomeFaFromFile(gf, getFaFile(opts));
  GenomeSam *gs = init_GenomeSam();
  loadGenomeSamFromFile(gs, getSamFile(opts));
  GenomeVcf_bplus *gv = genomeVcf_bplus_loadFile(getVcfFile(opts), 7, 6);

  const int cnt_thread = opt_threads(opts);
  const int64_t cnt_rec_sam = gsDataRecCnt(gs);
  pthread_t threads[cnt_thread];
  ThreadArgs args_thread[cnt_thread];

  // Assign arguments for threads
  for (int i = 0; i < cnt_thread; i++) {
    args_thread[i].id = i;
    args_thread[i].opts = opts;
    args_thread[i].gf = gf;
    args_thread[i].gs = gs;
    args_thread[i].gv = gv;
    args_thread[i].id_sam_start = (cnt_rec_sam / cnt_thread) * i;
  }
  for (int i = 0; i < cnt_thread - 1; i++) {
    args_thread[i].id_sam_end = args_thread[i + 1].id_sam_start - 1;
  }
  args_thread[cnt_thread - 1].id_sam_end = cnt_rec_sam - 1;

  // Create threads
  for (int i = 0; i < cnt_thread; i++) {
    if (pthread_create(&threads[i], NULL, integration_threads,
                       (void *)&args_thread[i]) != 0) {
      fprintf(
          stderr,
          "Error: failed to create thread to handle sam records from %" PRId64
          " to %" PRId64 "\n",
          args_thread->id_sam_start, args_thread->id_sam_end);
      exit(EXIT_FAILURE);
    }
  }

  // End threads
  for (int i = 0; i < cnt_thread; i++) {
    void *thread_ret = 0;
    pthread_join(threads[i], &thread_ret);
    printf("thread (%" PRId64 ") ended\n", (int64_t)thread_ret);
  }

  // Free structures

  destroy_GenomeFa(gf);
  destroy_GenomeSam(gs);
  destroy_GenomeVcf_bplus(gv);
  return;
}
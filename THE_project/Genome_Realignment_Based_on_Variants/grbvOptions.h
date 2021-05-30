/*
 * This header file must contain only simple defination and structures.
 * No subsequent process on defination and structures allowed.
 * @LastEditors: Atosh Dustosh
 */
#ifndef GRBVOPTIONS_H_INCLUDED
#define GRBVOPTIONS_H_INCLUDED

#pragma once

#include <stdio.h>
#include <stdlib.h>

/*
 * Operation types of Usage.
 */
#define NO_OPERATION 0

#define OPT_VERBOSE 1

/*
 * Set input and output fies.
 */
#define OPT_SET_OUTPUTFILE 101
#define OPT_SET_FAFILE 102
#define OPT_SET_FASTQFILE 103
#define OPT_SET_SAMFILE 104
#define OPT_SET_VCFFILE 105
#define OPT_SET_AUXFILE 106

static char *default_outputFile = "defaultOutput.txt";

/*
 * Set basic parameters.
 */

#define OPT_SET_SV_MIN_LEN 107
#define OPT_SET_SV_MAX_LEN 108

#define OPT_SET_MATCH 109
#define OPT_SET_MISMATCH 110
#define OPT_SET_GAPOPEN 111
#define OPT_SET_GAPEXTENSION 112

static const int default_sv_min_len = 51;
static const int default_sv_max_len = 300;

/*
 * Some simple operations against files.
 */
#define OPT_COUNTREC 201
#define OPT_FIRSTLINES 202
#define OPT_EXTRACTCHROM 203
#define OPT_STATISTICS_VCF 204

/*
 * GRBV operations.
 */
#define OPT_SELECTBADREADS 301
#define OPT_COMPARESAM 302

#define OPT_INTEGRATEVCFTOSAM 303
#define _OPT_INTEGRATION_SNPONLY 1
#define _OPT_INTEGRATION_SVONLY 2
#define _OPT_INTEGRATION_ALL 3

#define OPT_THREADS 401

#define OPT_KMERGENERATION 501

static const int default_kmerLength = 22;

#define OPTION_CONFLICT 1

typedef struct _define_Options {
  int ifOptConflict;
  int verbose;

  char *faFile;
  char *fastqFile;
  char *samFile;
  char *vcfFile;
  char *outputFile;
  char *auxFile;

  int sv_min_len;  // minimal length for a SV
  int sv_max_len;  // maximal length for a SV

  int match;
  int mismatch;
  int gapOpen;
  int gapExtension;

  int countRec;
  int firstLines;    // also store value of [firstline_number]
  int extractChrom;  // also store value of [chrom_idx]

  int selectBadReads;  // also store value of [MAPQ_threshold]
  int threads;         // also store value of [NUM_threads]

  int integration;  // also store selection of [integration_strategy]

  int kmerGeneration;  // also store value of [length_kmer]
} Options;

static inline void optCheck_conflict(Options *opts) {
  if (opts->ifOptConflict != OPTION_CONFLICT) {
    opts->ifOptConflict = OPTION_CONFLICT;
  } else {
    fprintf(stderr,
            "Warning: conflict options for this program. Please use multiple "
            "command lines to execute your tasks if needed, instead of "
            "running multiple operations in one command line. \n");
    exit(EXIT_FAILURE);
  }
}

/*
 * Methods for accessing data from a "Options *".
 */
static inline char *getFaFile(Options *opts) { return opts->faFile; }
static inline char *getFastqFile(Options *opts) { return opts->fastqFile; }
static inline char *getSamFile(Options *opts) { return opts->samFile; }
static inline char *getVcfFile(Options *opts) { return opts->vcfFile; }
static inline char *getOutputFile(Options *opts) { return opts->outputFile; }
static inline char *getAuxFile(Options *opts) { return opts->auxFile; }

static inline void setOutputFile(Options *opts, char *op_file) {
  opts->outputFile = op_file;
}

static inline int getSVminLen(Options *opts) { return opts->sv_min_len; }
static inline int getSVmaxLen(Options *opts) { return opts->sv_max_len; }
static inline int getMatch(Options *opts) { return opts->match; }
static inline int getMismatch(Options *opts) { return opts->mismatch; }
static inline int getGapopen(Options *opts) { return opts->gapOpen; }
static inline int getGapextension(Options *opts) { return opts->gapExtension; }

static inline int MAPQ_threshold(Options *opts) { return opts->selectBadReads; }

static inline int opt_threads(Options *opts) { return opts->threads; }

static inline int opt_integration_strategy(Options *opts) {
  return opts->integration;
}

static inline int opt_get_kmerLength(Options *opts) {
  return opts->kmerGeneration;
}
static inline void opt_set_kmerLength(int kmerLength, Options *opts) {
  opts->kmerGeneration = kmerLength;
}

typedef struct _define_FileList {
  char **paths;
  int count;
} FileList;

/**
 * @brief Get all files designated by command inputs.
 *
 * @param opts command inputs
 * @retval FileList* a list of designated files. You can access the value
 * "count" to get the size of it. The list must be freed mannually later with
 * destroyFileList().
 */
FileList *designatedFiles(Options *opts);

/**
 * @brief Destroy the FileList object
 */
void destroyFileList(FileList *fl);

#endif  // GRBVOPTIONS_H_INCLUDED

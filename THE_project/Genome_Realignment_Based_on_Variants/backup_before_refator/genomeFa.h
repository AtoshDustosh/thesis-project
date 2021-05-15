#ifndef RE_GENOMEFA_H_INCLUDED
#define RE_GENOMEFA_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "debug.h"
#include "genomeFaMacros.h"

/******************
 * Basic Structures
 ******************/

// This is actually a linked-list with an empty header
typedef struct _define_ChromFa {
  uint64_t *codedBases;  // binary coded bases using uint64_t array
  uint32_t length;       // length of chrom / number of bases (uncoded)
  char *info;            // info of chrom
  char *name;
  struct _define_ChromFa *next;
} ChromFa;

typedef struct _define_GenomeFa {
  uint16_t chromCnt;
  ChromFa *chroms;
} GenomeFa;

/**
 * @brief Initialize a ChromFa object and return the pointer to it. The
 * successfully returned object must be destroyed later using destroy_ChromFa()
 */
ChromFa *init_ChromFa();

/**
 * @brief Initialize a GenomeFa object and return the pointer to it. The
 * successfully returned object must be destroyed later using destroy_GenomeFa()
 */
GenomeFa *init_GenomeFa();

/**
 * @brief Destroy a ChromFa object. Should not be called by the user.
 */
void destroy_ChromFa(ChromFa *cf);

/**
 * @brief Destroy a GenomeFa object.
 */
void destroy_GenomeFa(GenomeFa *gf);

/************************************
 * Methods for manipulating GenomeFa
 ************************************/

/**
 * @brief  Get the chrom according to given info of that chromosome.
 * @retval pointer to the ChromFa object matching the give info; NULL if not
 * found
 */
ChromFa *getChromFromGenomeFabyInfo(const char *info, GenomeFa *gf);

/**
 * @brief  Similar to getChromFromGenomeFabyInfo, but this only use a substring
 * to recognize the chrom in GenomeFa, and thus is faster. This method is also
 * designed for synchronizing with sam-record-extraction. Recommend using this
 * to get a chrom.
 */
ChromFa *getChromFromGenomeFabyName(const char *name, GenomeFa *gf);

/**
 * @brief Get the chrom according to given index of that chromosome in the
 * GenomeFa object.
 *
 * @param idx 0-based index of the ChromFa object in the GenomeFa object. Note
 * that chrom 0 is the header chrom in GenomeFa object, and chrom 1, 2, ... are
 * the actual chroms containing bases. (so it is more like a 1-based index for
 * user)
 * @retval ChromFa* pointer to the ChromFa object matching the given info; NULL
 * if not found
 */
ChromFa *getChromFromGenomeFabyIndex(uint32_t idx, GenomeFa *gf);

/**
 * @brief  Get the base according to given position in the chromosome.
 * @param pos 1-based position of base
 * @retval Base type of the designated base
 */
Base getBase(ChromFa *cf, uint32_t pos);

/**
 * @brief  Get a specific base sequence from a ChromFa object.
 * @param  start: start position of the sequence in the chrom. 1-based position
 * and the base at the start position is included.
 * @param  end: end position of the sequence in the chrom. 1-based position and
 * the base at the end position is included.
 * @retval a string format of the base sequence. Must be freed later when not
 * needed anymore.
 */
char *getSeqFromChromFa(int64_t start, int64_t end, ChromFa *cf);

/**
 * @brief Add a ChromFa object into the GenomeFa object (linked to the end of
 * the linked-list with an empty header).
 */
void addChromToGenome(ChromFa *cf, GenomeFa *gf);

/**
 * @brief Load genome data into a GenomeFa object from designated file. Note
 * that the *.fa/*.fna file must match specifications.
 */
void loadGenomeFaFromFile(GenomeFa *gf, const char *filePath);

/**
 * @brief Write genome data into a designated file.
 */
void writeGenomeFaIntoFile(GenomeFa *gf, const char *filePath);

// ********************************
// Functions for Data Manipulating
// ********************************

/**
 * @brief  (helper function) Code the bpBuf array (a string with A,C,G,T) into a
 * uint64_t integer. Should not be called by the user. The calculation's
 * purpose: "ACGT" -> (0b) 001 010 011 100 000 000 ...... 000 0. Note that bases
 * start from the left and there may be some "0"s left out not filled. The bpBuf
 * string must contain no more than BP_PER_UINT64 chars (bases).
 */
static uint64_t codeBpBuf(char *bpBuf);

/**
 * @brief  Parse the info line from *.fa/*.fna file and fill in ChromFa.
 * @retval 0 for success; 1 for failure.
 */
static int parseFaInfo(ChromFa *cf, char *infoBuf);

/**********************************
 * Debugging Methods for GenomeFa
 **********************************/

void _testSet_genomeFa();

/**
 * @brief Print the GenomeFa object using terminal, but only print out the info
 * lines.
 */
void printGenomeFa_brief(GenomeFa *gf);

/**
 * @brief Print the GenomeFa object using terminal.
 */
void printGenomeFa(GenomeFa *gf);

/**
 * @brief Print the ChromFa object using terminal.
 */
void printChromFa(ChromFa *cf);

#endif
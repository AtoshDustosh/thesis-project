#ifndef GENOMEFA_H_INCLUDED
#define GENOMEFA_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "debug.h"
#include "genomeMacros.h"

// *****************
// Basic Structures
// *****************

// This is actually a doubly linked list with a header
// TODO note that if "cf->length" == 0, it may result in an ArrayOutOfBounder
// error. Planning to rewrite this macro.
#define getChromFaArrayLen(cf) ((cf->length - 1) / BP_PER_UINT64 + 1)

typedef struct _define_ChromFa {
  uint64_t *codedBases;          // binary coded bases using uint64_t array
  uint32_t length;               // length of chrom / number of bases (uncoded)
  char info[MAX_RECORD_LENGTH];  // info of chrom
  struct _define_ChromFa *prevChrom;
  struct _define_ChromFa *nextChrom;
} ChromFa;

typedef struct _define_GenomeFa {
  uint16_t chromNum;
  ChromFa *chroms;
} GenomeFa;

void _testSet_genomeFa();

/**
 * @brief Initialize a ChromFa object and return the pointer to it. Should not
 * be called by the user.
 */
static ChromFa *init_ChromFa();

/**
 * @brief Initialize a GenomeFa object and return the pointer to it.
 */
GenomeFa *init_GenomeFa();

/**
 * @brief  Get the chrom according to given info of that chromosome.
 * @retval pointer to the ChromFa object matching the give info; NULL if not
 * found
 */
static ChromFa *getChromFromGenome_info(char *info, GenomeFa *gf);

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
static ChromFa *getChromFromGenome_Idx(int idx, GenomeFa *gf);

/**
 * @brief  Get the base according to given position in the chromosome.
 * @param pos 1-based position of base
 * @retval Base type of the designated base
 */
Base getBase(ChromFa *cf, uint32_t pos);

/**
 * @brief Add a ChromFa object into the GenomeFa object (linked to the end of
 * the doubly-linked-list). Should not be called by the user.
 */
static void addChromToGenome(ChromFa *cf, GenomeFa *gf);

/**
 * @brief Remove a ChromFa object from the GenomeFa object. Should not be called
 * by the user.
 */
static void removeChromFromGenome(ChromFa *cf, GenomeFa *gf);

/**
 * @brief Destroy a ChromFa object. Should not be called by the user.
 */
static void destroy_ChromFa(ChromFa *cf);

/**
 * @brief Destroy a GenomeFa object.
 */
void destroy_GenomeFa(GenomeFa *gf);

/**
 * @brief Print the GenomeFa object using terminal.
 */
void printGenomeFa(GenomeFa *gf);

/**
 * @brief Print the GenomeFa object using terminal, but only print out the info
 * lines.
 */
void printGenomeFa_brief(GenomeFa *gf);

/**
 * @brief Print the ChromFa object using terminal.
 */
void printChromFa(ChromFa *cf);

// *****************
// Loading Functions
// *****************

/**
 * @brief  (helper function) Code the bpBuf array (a string with A,C,G,T) into a
 * uint64_t integer. Should not be called by the user. The calculation's
 * purpose: "ACGT" -> (0b) 001 010 011 100 000 000 ...... 000 0. Note that bases
 * start from the left and there may be some "0"s left out not filled. The bpBuf
 * string must contain no more than BP_PER_UINT64 chars (bases).
 */
static uint64_t codeBpBuf(char *bpBuf);

/**
 * @brief  (helper function) Receiving a new info line, pass it to this method.
 * It will create a new ChromFa object according to it, and the link it to the
 * end of GenomeFa->chroms. Should not be called by the user.
 * @retval 1 if infoBuf is valid and has successfully changed GenomeFa; 0
 * otherwise
 */
static int newInfoForGenomeFa(GenomeFa *gf, char *infoBuf);

/**
 * @brief Load genome data into a GenomeFa object from designated file. Note
 * that the *.fa/*.fna file must match specifications.
 */
void loadGenomeFaFromFile(GenomeFa *gf, FILE *fp);

/**
 * @brief Write genome data into a designated file.
 */
void writeGenomeFaIntoFile(GenomeFa *gf, FILE *fp);

#endif
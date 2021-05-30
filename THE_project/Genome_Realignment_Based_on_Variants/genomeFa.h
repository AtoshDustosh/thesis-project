#ifndef RE_GENOMEFA_H_INCLUDED
#define RE_GENOMEFA_H_INCLUDED

#pragma once

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "debug.h"
#include "genomeFaMacros.h"

/*********************************************************************
 *                     Structure Declarations
 ********************************************************************/

/**
 * @brief  A structure that keeps all information and bases of a chromosome.
 * @note   Although it is invalid to access this object' fields directly, you
 * can get this object from GenomeFa before extracting specific sequences on the
 * chromosome. And thus you can use the object to extract sequences without
 * searching for it again.
 */
typedef struct ChromFa ChromFa;

typedef struct GenomeFa GenomeFa;

/**
 * @brief Initialize a GenomeFa object and return the pointer to it. The
 * successfully returned object must be destroyed later using destroy_GenomeFa()
 */
GenomeFa *init_GenomeFa();

/**
 * @brief Destroy a GenomeFa object.
 */
void destroy_GenomeFa(GenomeFa *gf);

/*********************************************************************
 *                           Data Accessor
 ********************************************************************/

/**
 * @brief  Get length of the chromosome.
 */
extern uint32_t chromFa_length(ChromFa *cf);

/**
 * @brief  Get name of the chromosome.
 */
extern const char *chromFa_name(ChromFa *cf);

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
 * @param idx 1-based index of the ChromFa object in the GenomeFa object. Note
 * that chrom 0 is the header chrom in GenomeFa object, and chrom 1, 2, ... are
 * the actual chroms containing bases.
 * @retval ChromFa* pointer to the ChromFa object matching the given info; NULL
 * if not found
 */
ChromFa *getChromFromGenomeFabyIndex(uint32_t idx, GenomeFa *gf);

/**
 * @brief  Get the base according to given position in the chromosome.
 * @param  pos 1-based position of base
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
 * @brief  Transit local position on a chromsome into the absolute position on the whole genome. For example, the current absolute position equals to the offset of current chromosome plus the ending position of its last chromosome.
 * @param  id_chrom: 1-based index for chromsomes.
 * @param  pos: 1-based position on current chromosome.
 * @retval absolute position on the whole genome
 */
int32_t genomeFa_absolutePos(uint32_t id_chrom, int32_t pos, GenomeFa *gf);

/*********************************************************************
 *                      Data Loading and Writing
 ********************************************************************/

/**
 * @brief  Load genome data into a GenomeFa object from designated file.
 */
GenomeFa *genomeFa_loadFile(char *filePath);

/**
 * @brief Write genome data into a designated file.
 */
void writeGenomeFaIntoFile(GenomeFa *gf, const char *filePath);

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
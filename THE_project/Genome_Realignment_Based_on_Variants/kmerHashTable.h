#ifndef KMERHASHTABLE_H_INCLUDED
#define KMERHASHTABLE_H_INCLUDED

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/*********************************************************************
 *                         Structures and Accessors
 ********************************************************************/

typedef struct KmerHashCell KmerHashCell;

typedef struct KmerHashTable KmerHashTable;

typedef struct Iterator_Kmerhash Iterator_Kmerhash;

extern char *kmerHashCell_string(KmerHashCell *cell);

extern int32_t kmerHashCell_pos(KmerHashCell *cell);

extern KmerHashCell *kmerHashCell_next(KmerHashCell *cell);

/*********************************************************************
 *                            Basic Functions
 ********************************************************************/

Iterator_Kmerhash* init_iterator_Kmerhash(KmerHashTable *table);

/**
 * @brief  Get a copy of temporary string (kmer). 
 */
char *iterator_Kermhash_string(Iterator_Kmerhash *it);

/**
 * @brief  Get temporary position of the kmer
 */
int32_t iterator_Kmerhash_pos(Iterator_Kmerhash *it);

/**
 * @brief  Do iteration with the iterator. 
 * @retval true if there exists another record; false otherwise.
 */
bool iterator_Kmerhash_next(Iterator_Kmerhash *it);

void destroy_iterator_Kmerhash(Iterator_Kmerhash *it);

KmerHashTable *init_kmerHashTable(int32_t tableSize, int32_t kmerLength);

void destroy_kmerHashTable(KmerHashTable *table);

/**
 * @brief  Add a data unit (kmerString, kmerPos) into the hash table
 * @param  *kmer: string of the kmer (char array ended with '\0')
 * @param  pos: 1-based position of the kmer
 * @retval true if succeeded; false if there exists identical kmer
 */
bool kmerHashTable_add(char *kmer, int32_t pos, KmerHashTable *table);

/**
 * @brief  Check the existence of a data unit (kmerString, kmerPos).
 * @param  *kmer: string of the kmer (char array ended with '\0')
 * @param  pos: 1-based position of the kmer
 * @retval true if exists the same kmer; false otherwise
 */
bool kmerHashTable_exist(char *kmer, int32_t pos, KmerHashTable *table);

/*********************************************************************
 *                            Debug Functions
 ********************************************************************/

/**
 * @brief  Print every hash cell and all data within it.
 */
void print_kmerHashTable(KmerHashTable *table);

/**
 * @brief  Collect and print statistical data of the table.
 */
void statistics_kmerHashTable(KmerHashTable *table);

int _testSet_hashTable();

#endif
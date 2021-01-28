#ifndef HASHTABLE_H_INCLUDED
#define HASHTABLE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

typedef struct _define_HashCell {
    uint64_t data;
    int flag;   // used to count how many hash cells are after it
    struct _define_HashCell* nextCell;
} HashCell;

typedef struct _define_HashTable {
    HashCell** hashList;
    uint64_t tableSize;
} HashTable;


/**
 * A collection of test in this header file.
 */
void _HashTableTestSet();


/*
 * Working functions.
 */


/**
 * Get the first hash cell of the hash index of the string if the string is to be added into the
 * hash table.
 *
 * @param str string to be searched in hash table
 * @param hashTable hash table
 * @param tableSize size of the hash table
 * @return first hash cell of the hash index of the input string
 */
HashCell* searchHashIndexOfString(char* str, HashTable* hashTable, uint64_t tableSize);

/**
 * Check performance of the hash table.
 *
 * @param hashTable hash table
 * @param tableSize size of the table
 */
void checkHashTablePerformance(HashTable* hashTable, uint64_t tableSize);

/**
 * Display the hash table.
 *
 * @param hashTable hash table
 * @param tableSize size of the table
 */
void displayHashTable(HashTable* hashTable, uint64_t tableSize);

/**
 * Add a hash string with specific data as a hash cell into hash table.
 *
 * @param str string used for calculating hash index
 * @param data data that the string represents
 * @param hashTable hash table
 * @param tableSize size of the hash table
 */
void addHashCell(char* str, uint64_t data, HashTable* hashTable, uint64_t tableSize);

/**
 * Initialize a hash table of specific size with all hash cells being NULL.
 *
 * @param tableSize size of the hash table
 * @return hash table
 */
HashTable* initHashTable(uint64_t tableSize);

/**
 * Calculating the hash table index of a string.
 *
 * @param str input string
 * @param tableSize size of the hash table
 * @return hash table index
 */
uint64_t hashIndex(char *str, uint64_t tableSize);

#endif // HASHTABLE_H_INCLUDED

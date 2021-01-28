#ifndef FILEOPERATION_H_INCLUDED
#define FILEOPERATION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "AuxiliaryDataType.h"
#include "Read.h"

/**
 * Get the size of data (DNA sequence) of a ?.fna file.
 *
 * @param filePath file path
 * @return data size of ?.fna file
 */
uint64_t fnaDataSize(char* filePath);

/**
 * Load ?.fna data file into memory - stored in an uint64_t[] array.
 *
 * @param filePath file path
 * @param dataLength length of ?.fna data
 * @param / @return hexCodedDNA hex-coded DNA sequence
 * @param / @return file header of ?.fna data file
 */
void loadFnaData(char* filePath, uint64_t dataLength, uint64_t* hexCodedDNA, char* fnaFileHeader);

/**
 * Open a ?.fastq file and load a read into memory.
 * If fpointer is not NULL, directly load a read into memory.
 *
 * @param filePath file path
 * @param fpointer pointer to a file pointer - points to start position in file
 *      of next read when the method finishes
 * @param read Read data type used for storing a read
 * @param 1 if there are more reads to be processed in the ?.fastq file; 0 otherwise
 */
uint64_t loadOneReadFromFile(char* filePath, FILE** fpointer, Read* read);


#endif // FILEOPERATION_H_INCLUDED

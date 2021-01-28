#ifndef SNAP_H_INCLUDED
#define SNAP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "HashTable.h"
#include "Read.h"
#include "AuxiliaryDataType.h"


/*
 * Results of SNAP alignment.
 */
#define MULTIPLE_HITS 2
#define SINGLE_HIT 1
#define NOT_FOUND 0



/**
 * A struct type for SNAP.
 */
typedef struct _define_SNAP {
    uint64_t* hexCodedRefDNA; // hex-coded reference DNA
    uint64_t DNAlength;       // length of DNA
    uint64_t seedLength;      // length of a seed
    HashTable* hashTable;     // hash table constructed out of reference DNA
    Read* read;              // a read (segment) - SNAP only process one read per time
} SNAP;






/**
 * A collection of test in this header file.
 */
void _SNAPTestSet();




/*
 * Working functions.
 */


/**
 * Construct the basic information of a SNAP structure.
 *
 * @param hexCodedRefDNA hex-coded reference DNA
 * @param DNAlength length of DNA
 * @param seedLength length of seed
 * @return an SNAP structure
 */
SNAP* constructSNAP(uint64_t* hexCodedRefDNA, uint64_t DNAlength, uint64_t seedLength);

/**
 * Load a read into a SNAP structure.
 *
 * @param read a read
 * @param snap an SNAP structure
 */
void loadOneReadIntoSNAP(Read* read, SNAP* snap);

/**
 * Create a new hex-coded string buffer.
 * Extract a fragment of specific length from reference sequence and put it into the hex-coded string
 * buffer.
 *
 * @param refSeq hex-coded reference sequence
 * @param refLength length of reference sequence
 * @param fragmentLength length of fragment
 * @param refOffset offset of the fragment on ref sequence
 * @return hex-coded string buffer of the fragment
 */
HexCodedStringBuffer* extractHexCodedFragmentFromRef(uint64_t* refSeq, uint64_t refLength,
        uint64_t fragmentLength, uint64_t refOffset);


/**
 * Align a read using SNAP with specific parameters.
 *
 * @param snap an SNAP structure that is already fully constructed
 * @param seedLength length of seed (must be <= 32)
 * @param EDmax maximum edit-distance
 * @param hitMax maximum hit count per seed
 * @param confidenceThreshold confidence threshold - limit the difference of the best alignment
 *      and the second-best alignment
 * @return single hit, multiple hits or not found
 */
uint64_t alignOneReadUsingSNAP(SNAP* snap, uint64_t seedLength, uint64_t EDmax,
                             uint64_t hitMax, uint64_t confidenceThreshold);









#endif // SNAP_H_INCLUDED

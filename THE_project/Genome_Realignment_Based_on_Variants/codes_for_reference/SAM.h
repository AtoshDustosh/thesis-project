#ifndef SAM_H_INCLUDED
#define SAM_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <inttypes.h>

#include "Read.h"

typedef struct _define_SAM {
    char* filePath;
    FILE* fpointer;
    Read* lastRead;
    Read* newRead;
} SAM;

/**
 * Construct a SAM structure.
 *
 * @param char* filePath path of ?.sam file
 * @param refName name of reference DNA
 * @param refLength length of reference DNA
 * @return a new SAM structure
 */
SAM* constructSAM(char* filePath, char* refName, uint64_t refLength);


/**
 * Print an aligned read using a SAM structure.
 *
 * @param read a new aligned read
 * @param sam a SAM structure
 */
void fprintReadUsingSAM(Read* read, SAM* sam);



#endif // SAM_H_INCLUDED

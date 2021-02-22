#ifndef READ_H_INCLUDED
#define READ_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

/**
 * A type used for storing the alignment results of a read.
 *
 */
typedef struct _define_Read {
    char* QNAME; // name of read
    int FLAG;           // bit-flags
    char* RNAME; // name of reference
    int POS;            // 1-based position = offset + 1
    float MAPQ;         // (pass)
    char* CIGAR; // CIGAR (M, D, I)
    char* RNEXT; // name of reference of mate read alignment
    int PNEXT;          // offset on ref of mate read
    int TLEN;           // length of template
    char* SEQ;   // sequence of read (A, C, G, T)
    char* QUAL;  // quality of read (Q = -10logE)
} Read;


/**
 * Initialize a struct Read.
 *
 * @param read a Read struct type
 */
void initRead(Read *read);


/**
 * Print all fields of a read.
 *
 * @param read a Read struct type
 */
void printRead(Read* read);

/**
 * Get a copy of Read.
 *
 * @param read a Read struct type
 * @return a copy of read
 */
Read* copyRead(Read* read);

/**
 * Clear a read and free its memory.
 *
 * @param read a Read struct type
 */
void clearRead(Read* read);

#endif // READ_H_INCLUDED

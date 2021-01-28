#ifndef LANDAUVISHKIN_H_INCLUDED
#define LANDAUVISHKIN_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "AuxiliaryDataType.h"

/**
 * A collection of test in this header file.
 */
void _EditDistanceTestSet();

/*
 * Working functions.
 */

/**
 * Calculate edit distance between 2 strings.
 *
 * @param patternStrBuf string buffer of pattern sequence
 * @param refStrBuf string buffer of reference sequence
 * @param EDmax maximum edit distance that can be accepted
 * @param / @return CIGARbuffer buffer for CIGAR string
 * @param maxBufLen maximum length of buffer for CIGAR string
 * @return best edit-distance if it's within limit; INITEDVALUE if out of limit
 */
uint64_t calculateEditDistance(StringBuffer* patternStrBuf, StringBuffer* refStrBuf, uint64_t EDmax,
                               char* CIGARbuffer, uint64_t maxBufLen);

#endif // LANDAUVISHKIN_H_INCLUDED

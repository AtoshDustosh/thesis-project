#ifndef AUXILIARYFUNCTION_H_INCLUDED
#define AUXILIARYFUNCTION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "AuxiliaryDataType.h"

/**
 * A collection of test in this header file.
 */
void _AuxiliaryFunctionTestSet();





/*
 * Working functions.
 */


/**
 * Copy a string and return a pointer to the copy.
 *
 * @param str the string to be copied
 * @return a copy of the string
 */
char* copyString(char* str);

/**
 * Reverse a string.
 *
 * @param / @return str string to be reversed
 */
void reverseString(char* str);

/**
 * Transform a hex-coded string buffer to a string buffer.
 *
 * @param hexCodedStrBuf hex-coded string buffer
 * @param / @return strBuf string buffer
 * @param charNumPerHex #(chars) per hexadecimal number
 */
 StringBuffer* transHexCodedStringBufferToStringBuffer(HexCodedStringBuffer* hexCodedStrBuf,
                                                       uint64_t charNumPerHex);

/**
 * Extract the hexadecimal bits of a char compressed in a 64-bit hexInt.
 *
 * @param offset offset of in the 64-bit hexInt (equals to index of the char)
 * @param hexInt a 64-bit hexadecimal integer
 * @param charNumPerHex #(chars) per hexadecimal number
 * @return extracted hexadecimal number to the corresponding index
 */
uint64_t extractCharBitFromHexInt(uint64_t offset, uint64_t hexInt, uint64_t charNumPerHex);

/**
 * Transform a string buffer to a hex-coded string buffer.
 *
 * @param strBuf string buffer
 * @param / @return hexCodedStrBuf hex-coded string buffer
 * @param charNumPerHex #(chars) per hexadecimal number
 */
HexCodedStringBuffer* transStringBufferToHexCodedStringBuffer(StringBuffer* strBuf,
                                                              uint64_t charNumPerHex);

/**
 * Get the lower case of a character.
 *
 * @param ch an English letter
 * @return lower case of ch if ch is in upper case, i.e. 'A' -> 'a';
 *      ch otherwise, i.e. 'a' -> 'a', '*' -> '*'
 */
char lowerCase(char ch);

/**
 * Get the upper case of a character.
 *
 * @param ch an English letter
 * @return upper case of ch if ch is in lower case, i.e. 'a' -> 'A';
 *      ch otherwise, i.e. 'A' -> 'A', '*' -> '*'
 */
char UpperCase(char ch);

/**
 * Transform hexadecimal numbers into characters (a, c, g, t)
 *
 * @param hexValue hexadecimal number
 * @return 'a' if 0x0; 'c' if 0x1;
 *         'g' if 0x2; 't' if 0x3;
 *         '*' otherwise.
 */
char hexToChar(uint64_t hex);

/**
 * Transform characters (A, C, G, T) into hexadecimal numbers.
 *
 * @param ch character
 * @return 0x0 if ch == 'A' || ch == 'a';
 *         0x1 if ch == 'C' || ch == 'c';
 *         0x2 if ch == 'G' || ch == 'g';
 *         0x3 if ch == 'T' || ch == 't';
 *         0x0 otherwise.
 */
uint64_t charToHex(char ch);

/**
 * Get the hex-code of the inverse of a base.
 *
 * @param base base
 * @return hex-code of the inverse of a base
 */
uint64_t getInverseBaseHex(uint64_t base);

/**
 * Get the min value between 2 (uint64_t type) values.
 *
 * @param value1 value1
 * @param value2 value2
 */
uint64_t min_uint64_t(uint64_t value1, uint64_t value2);

/**
 * Get the max value between 2 (uint64_t type) values.
 *
 * @param value1 value1
 * @param value2 value2
 */
uint64_t max_uint64_t(uint64_t value1, uint64_t value2);

#endif // AUXILIARYFUNCTION_H_INCLUDED

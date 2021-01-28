#ifndef AUXILIARYDATATYPE_H_INCLUDED
#define AUXILIARYDATATYPE_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>


#define HEX_CODED_STRINGBUFFER_SAME 0
#define HEX_CODED_STRINGBUFFER_DIFFERNET 1

/**
 * A type used for storing a string.
 *
 * \note need to apply cpu memory for its buffer before usage
 */
typedef struct _define_StringBuffer {
    char* buffer;    // buffer used to store a string
    uint64_t length; // length of the string stored in this StringBuffer
} StringBuffer;

/**
 * A type used for storing a hexadecimal-coded string.
 *
 * \note a string is coded into one or multiple 64-bit unsigned hexadecimals and stored
 * in this HexadecimalStringBuffer. This type will store all hexadecimals, length of the
 * hex-array, and length of the string.
 *
 * \note need to apply cpu memory for its buffer before usage
 */
typedef struct _define_HexCodedStringBuffer {
    uint64_t* hexArray; // an array of 64-bit hexadecimal numbers
    uint64_t arrayLength;   // length of the hexArray
    uint64_t strLength; // length of the string i.e. #(characters or letters)
} HexCodedStringBuffer;







/**
 * A collection of test in this header file.
 */
void _AuxiliaryDataTypeTestSet();



/*
 * Working functions.
 */

/**
 * Construct a hex-coded string buffer with input (hexArray, arrayLength, strLength) by copying
 * hexArray.
 *
 * @param hexCodedStrBuf hex-coded string buffer to be constructed - cannot be NULL pointer
 * @param hexArray hex-coded array
 * @param arrayLength length of array
 * @param strLength length of string
 */
HexCodedStringBuffer* constructHexCodedStringBuffer(uint64_t* hexArray, uint64_t arrayLength,
                                                    uint64_t strLength);

/**
 * Construct a string buffer with input (buffer, length) by copying buffer.
 *
 * @param strBuf string buffer to be constructed - cannot be NULL pointer
 * @param buffer string
 * @param length string length
 */
StringBuffer* constructStringBuffer(char* buffer, uint64_t length);

/**
 * Compare 2 hex-coded string buffer.
 *
 * @param hexCodedStrBuf1 hex-coded string buffer 1
 * @param hexCodedStrBuf2 hex-coded string buffer 2
 * @return HEX_CODED_STRINGBUFFER_SAME if hex-buffer-1 == hex-buffer-2;
 *          HEX_CODED_STRINGBUFFER_DIFFERENT otherwise
 */
uint64_t compareHexCodedStringBuffer(HexCodedStringBuffer* hexCodedStrBuf1,
                                     HexCodedStringBuffer* hexCodedStrBuf2);


/**
 * Print a string buffer.
 *
 * @param strBuf string buffer - cannot be NULL pointer
 */
void printStringBuffer(StringBuffer* strBuf);

/**
 * Print a hex-coded string buffer.
 *
 * @param hexCodedStrBuf hex-coded string buffer - cannot be NULL pointer
 */
void printHexCodedStringBuffer(HexCodedStringBuffer* hexCodedStrBuf);

/**
 * Clear the memory occupation of a string buffer.
 *
 * @param strBuf a string buffer
 */
void clearStringBuffer(StringBuffer* strBuf);

/**
 * Clear the memory occupation of a hex-coded string buffer.
 *
 * @param hexCodedStrBuf a hex-coded string buffer
 */
void clearHexCodedStringBuffer(HexCodedStringBuffer* hexCodedStrBuf);

#endif // AUXILIARYDATATYPE_H_INCLUDED

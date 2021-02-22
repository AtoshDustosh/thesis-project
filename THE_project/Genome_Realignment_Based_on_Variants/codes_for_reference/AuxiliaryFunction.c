#include "AuxiliaryFunction.h"

#include <math.h>
#include <string.h>

#include "MyArgs.h"

#define HEX_FOR_LETTER_A 0x0
#define HEX_FOR_LETTER_C 0x1
#define HEX_FOR_LETTER_G 0x2
#define HEX_FOR_LETTER_T 0x3


static void _getInverseBaseHexTest();
static void _charToHexTest();
static void _hexToCharTest();
static void _lowerCaseTest();
static void _transStringBufferToHexCodedStringBufferTest();
static void _extractCharBitsFromHexIntTest();
static void _transHexCodedStringBufferToStringBufferTest();


void _AuxiliaryFunctionTestSet() {
    _getInverseBaseHexTest();
    _charToHexTest();
    _hexToCharTest();
    _lowerCaseTest();
    _transStringBufferToHexCodedStringBufferTest();
    _extractCharBitsFromHexIntTest();
    _transHexCodedStringBufferToStringBufferTest();
}


/*
 * Tests for working functions.
 */





/**
 * Test function getInverseBaseHex.
 */
static void _getInverseBaseHexTest() {
    uint64_t hexA = HEX_FOR_LETTER_A;
    uint64_t hexC = HEX_FOR_LETTER_C;
    uint64_t hexG = HEX_FOR_LETTER_G;
    uint64_t hexT = HEX_FOR_LETTER_T;

    printf("\n**************** _getInverseBaseHexTest ****************\n");
    printf("%c, 0x%"PRIx64" -> %c, 0x%"PRIx64"\n", hexToChar(hexA), hexA,
           hexToChar(getInverseBaseHex(hexA)), getInverseBaseHex(hexA));
    printf("%c, 0x%"PRIx64" -> %c, 0x%"PRIx64"\n", hexToChar(hexC), hexC,
           hexToChar(getInverseBaseHex(hexC)), getInverseBaseHex(hexC));
    printf("%c, 0x%"PRIx64" -> %c, 0x%"PRIx64"\n", hexToChar(hexG), hexG,
           hexToChar(getInverseBaseHex(hexG)), getInverseBaseHex(hexG));
    printf("%c, 0x%"PRIx64" -> %c, 0x%"PRIx64"\n", hexToChar(hexT), hexT,
           hexToChar(getInverseBaseHex(hexT)), getInverseBaseHex(hexT));
}

/**
 * Test function charToHex.
 */
static void _charToHexTest() {
    printf("\n**************** _charToHexTest ****************\n");
    char a = 'a';
    char A = 'A';
    char c = 'c';
    char C = 'C';
    char g = 'g';
    char G = 'G';
    char t = 't';
    char T = 'T';
    printf("%c, 0x%"PRIx64" | %c, 0x%"PRIx64"\n", a, charToHex(a), A, charToHex(A));
    printf("%c, 0x%"PRIx64" | %c, 0x%"PRIx64"\n", c, charToHex(c), C, charToHex(C));
    printf("%c, 0x%"PRIx64" | %c, 0x%"PRIx64"\n", g, charToHex(g), G, charToHex(G));
    printf("%c, 0x%"PRIx64" | %c, 0x%"PRIx64"\n", t, charToHex(t), T, charToHex(T));
    printf("%c, 0x%"PRIx64"\n", '$', charToHex('$'));
    printf("%c, 0x%"PRIx64"\n", '2', charToHex('2'));
}

/**
 * Test function hexToChar.
 */
static void _hexToCharTest() {
    printf("\n**************** _hexToCharTest ****************\n");
    uint64_t hexA = HEX_FOR_LETTER_A;
    uint64_t hexC = HEX_FOR_LETTER_C;
    uint64_t hexG = HEX_FOR_LETTER_G;
    uint64_t hexT = HEX_FOR_LETTER_T;
    uint64_t hexWTF = 0xa;
    printf("0x%"PRIx64" -> %c\n", hexA, hexToChar(hexA));
    printf("0x%"PRIx64" -> %c\n", hexC, hexToChar(hexC));
    printf("0x%"PRIx64" -> %c\n", hexG, hexToChar(hexG));
    printf("0x%"PRIx64" -> %c\n", hexT, hexToChar(hexT));
    printf("0x%"PRIx64" -> %c\n", hexWTF, hexToChar(hexWTF));
}

/**
 * Test function lowerCase.
 */
static void _lowerCaseTest() {
    printf("\n**************** _lowerCaseTest ****************\n");
    printf("lowerCase(%c) -> %c\n", 'A', lowerCase('A'));
    printf("lowerCase(%c) -> %c\n", 'a', lowerCase('a'));
    printf("lowerCase(%c) -> %c\n", '#', lowerCase('#'));
}

/**
 * Test function transBufToHexInt.
 */
static void _transStringBufferToHexCodedStringBufferTest() {
    printf("\n*************** _transStringBufferToHexCodedStringBufferTest ***************\n");
    const uint64_t charNumPerHex = CHAR_NUM_PER_HEX;
    StringBuffer* strBuf = NULL;
    HexCodedStringBuffer* hexCodedStrBuf = NULL;

    printf("... initialization completed \n");


// test string length within charNumPerHex
    // 01010101 00100010 01110111 01010101 00000000 10101010 11111111 00000101
    // 0x5522775500AAFF05
    strBuf = constructStringBuffer("ccccagagctctccccaaaaggggttttaacc", 32);
    hexCodedStrBuf = transStringBufferToHexCodedStringBuffer(strBuf, charNumPerHex);
    printStringBuffer(strBuf);
    printHexCodedStringBuffer(hexCodedStrBuf);
    printf("expected: 0x5522775500aaff05\n");
    clearStringBuffer(strBuf);
    clearHexCodedStringBuffer(hexCodedStrBuf);

    // 01010101 00100010 01110111 01010101 00000000 00000000 00000000 00000000
    // 0x5522775500AAFF00
    strBuf = constructStringBuffer("ccccagagctctcccc", 16);
    hexCodedStrBuf = transStringBufferToHexCodedStringBuffer(strBuf, charNumPerHex);
    printStringBuffer(strBuf);
    printHexCodedStringBuffer(hexCodedStrBuf);
    printf("expected: 0x5522775500000000\n");
    clearStringBuffer(strBuf);
    clearHexCodedStringBuffer(hexCodedStrBuf);

// test string length over charNumPerHex
    // 01010101 00100010 01110111 01010101 00000000 10101010 11111111 00000000 ~
    // ~ 11001100 10001000 11001100 01000100
    // 0x5522775500AAFF00 0xCC88CC4400000000
    strBuf = constructStringBuffer("ccccagagctctccccaaaaggggttttaaaatatagagatatacaca", 48);
    hexCodedStrBuf = transStringBufferToHexCodedStringBuffer(strBuf, charNumPerHex);
    printStringBuffer(strBuf);
    printHexCodedStringBuffer(hexCodedStrBuf);
    printf("expected: 0x5522775500aaff00, 0xcc88cc4400000000\n");
    clearStringBuffer(strBuf);
    clearHexCodedStringBuffer(hexCodedStrBuf);
}

/**
 * Test function extractCharBitsFromHexInt.
 */
static void _extractCharBitsFromHexIntTest() {
    printf("\n**************** _extractCharBitsFromHexIntTest ****************\n");
    uint64_t hexInt = 0x27fd3de1e41a90ce;
    char* charSequence = "agcttttcattctgactgcaacgggcaatatg";
    uint64_t i = 0;
    uint64_t length = CHAR_NUM_PER_HEX;

    printf("0x%16"PRIx64"\n", hexInt);
    printf("%s\n", charSequence);
    for(i = 0; i < length; i++) {
        int offset = i % CHAR_NUM_PER_HEX;
        uint64_t extractedHex = extractCharBitFromHexInt(offset, hexInt, CHAR_NUM_PER_HEX);
        printf("%c-0x%"PRIx64" ", charSequence[i], extractedHex);
        if((i + 1) % 8 == 0) {
            printf("\n");
        }
    }
}

/**
 * Test function transHexCodedStringBufferToStringBuffer.
 */
static void _transHexCodedStringBufferToStringBufferTest() {
    printf("\n*************** _transHexCodedStringBufferToStringBufferTest ***************\n");
    const uint64_t charNumPerHex = CHAR_NUM_PER_HEX;

    /*
     * Test single-hexInt string buffer transformation.
     */
    StringBuffer* strBuf = NULL;
    uint64_t arrayLength = 1;
    uint64_t strLength = 32;
    uint64_t* hexArray1 = (uint64_t*)malloc(sizeof(uint64_t) * arrayLength);
    if(hexArray1 == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    *(hexArray1) = 0x27fd3de1e41a90ce;
    HexCodedStringBuffer* hexCodedStrBuf =
        constructHexCodedStringBuffer(hexArray1, arrayLength, strLength);

    printHexCodedStringBuffer(hexCodedStrBuf);
    strBuf = transHexCodedStringBufferToStringBuffer(hexCodedStrBuf, charNumPerHex);
    printStringBuffer(strBuf);
    printf("expected: agcttttcattctgactgcaacgggcaatatg -> strcmp:%d\n",
           strcmp(strBuf->buffer, "agcttttcattctgactgcaacgggcaatatg"));
    clearStringBuffer(strBuf);
    clearHexCodedStringBuffer(hexCodedStrBuf);

    /*
     * Test multiple-hexInt string buffer transformation.
     */
    /** < \note if you use "uint64 array[2]" to create an array and pass it to hex-coded string
        buffer. Freeing the array in hex-coded string buffer will result in bugs. */
    arrayLength = 2;
    strLength = 57;
    uint64_t* hexArray2 = (uint64_t*)malloc(sizeof(uint64_t) * arrayLength);
    if(hexArray2 == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    *(hexArray2) = 0x27fd3de1e41a90ce;
    *(hexArray2 + 1) = 0xffff0000ffff0000;
    hexCodedStrBuf = constructHexCodedStringBuffer(hexArray2, arrayLength, strLength);

    printHexCodedStringBuffer(hexCodedStrBuf);
    strBuf = transHexCodedStringBufferToStringBuffer(hexCodedStrBuf, charNumPerHex);
    printStringBuffer(strBuf);
    printf("expected: agcttttcattctgactgcaacgggcaatatgttttttttaaaaaaaatttttttta -> strcmp:%d\n",
           strcmp(strBuf->buffer, "agcttttcattctgactgcaacgggcaatatgttttttttaaaaaaaatttttttta"));
    clearStringBuffer(strBuf);
    clearHexCodedStringBuffer(hexCodedStrBuf);
}

















/*
 * Working functions.
 */

char* copyString(char* str) {
    char* copiedStr = (char*)malloc(sizeof(char) * strlen(str));
    if(copiedStr == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    uint64_t i = 0;
    while(*(str + i) != '\0'){
        copiedStr[i] = str[i];
        i++;
    }
    copiedStr[i] = '\0';
    return copiedStr;
}



void reverseString(char* str) {
    uint64_t strLength = (uint64_t)strlen(str);
    char* temp = (char*)malloc(sizeof(char) * (strLength + 1));
    if(temp == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    for(uint64_t i = 0; i < strLength; i++) {
        temp[i] = str[i];
    }
    temp[strLength] = '\0';

    for(uint64_t i = 0; i < strLength; i++) {
        str[i] = temp[strLength - 1 - i];
    }
    free(temp);
}

StringBuffer* transHexCodedStringBufferToStringBuffer(HexCodedStringBuffer* hexCodedStrBuf,
        uint64_t charNumPerHex) {
    StringBuffer* strBuf = (StringBuffer*)malloc(sizeof(StringBuffer));
    if(strBuf == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    strBuf->buffer = (char*)malloc(sizeof(char) * (hexCodedStrBuf->strLength + 1));
    if(strBuf->buffer == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    strBuf->length = hexCodedStrBuf->strLength;
    uint64_t* hexArray = hexCodedStrBuf->hexArray;

    uint64_t i = 0;

    for(i = 0; i < hexCodedStrBuf->strLength; i++) {
        uint64_t hexArrayIndex = i / charNumPerHex;
        uint64_t hexInt = *(hexArray + hexArrayIndex);
        uint64_t charHex = extractCharBitFromHexInt(i, hexInt, charNumPerHex);
        strBuf->buffer[i] = hexToChar(charHex);
//        printf("charHex: 0x%"PRIx64", char: %c\n", charHex, strBuf->buffer[i]);
    }
    // malloc will not append a '\0' to the end of a char*.
    strBuf->buffer[i] = '\0';
    return strBuf;
}

uint64_t extractCharBitFromHexInt(uint64_t offset, uint64_t hexInt, uint64_t charNumPerHex) {
    uint64_t bitInterval = sizeof(uint64_t) * 8 / charNumPerHex;
    uint64_t bitShiftLeft = offset * bitInterval;
    uint64_t bitShiftRight = bitInterval * (charNumPerHex - 1);

    return (hexInt << bitShiftLeft) >> bitShiftRight;
}

HexCodedStringBuffer* transStringBufferToHexCodedStringBuffer(StringBuffer* strBuf,
        uint64_t charNumPerHex) {
    HexCodedStringBuffer* hexCodedStrBuf =
        (HexCodedStringBuffer*)malloc(sizeof(HexCodedStringBuffer));
    if(hexCodedStrBuf == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    // initialize hex-coded string buffer
    if(strBuf->length % charNumPerHex == 0) {
        hexCodedStrBuf->arrayLength = strBuf->length / charNumPerHex;
    } else {
        hexCodedStrBuf->arrayLength = strBuf->length / charNumPerHex + 1;
    }

    hexCodedStrBuf->hexArray =
        (uint64_t*)malloc(sizeof(uint64_t) * hexCodedStrBuf->arrayLength);
    if(hexCodedStrBuf->hexArray == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    hexCodedStrBuf->strLength = strBuf->length;

    uint64_t bitInterval = sizeof(uint64_t) * 8 / charNumPerHex;
    uint64_t str_index = 0;
    uint64_t array_i = 0;
    char bufChar = strBuf->buffer[str_index++];
    uint64_t hexInt = 0x0;
    while(bufChar != '\0') {
        uint64_t hexBit = charToHex(bufChar);
        uint64_t bitShift = (charNumPerHex - (str_index % charNumPerHex)) % charNumPerHex;
        hexInt = hexInt | (hexBit << (bitShift * bitInterval));
//        printf("%c\t%"PRIu64"\t%#"PRIx64"\n", bufChar, bitShift, hexInt);
        bufChar = strBuf->buffer[str_index++];
        // put hexInt into hex-coded string buffer when it's fully coded
        if(bitShift == 0 || bufChar == '\0') {
            hexCodedStrBuf->hexArray[array_i++] = hexInt;
//            printf("########hexInt: %#"PRIx64"\n", hexInt);
            hexInt = 0x0;
        }
    }
    return hexCodedStrBuf;
}

char lowerCase(char ch) {
    if(ch >= 'A' && ch <= 'Z') {
        return ch - ('A' - 'a');
    } else {
        return ch;
    }
}

char UpperCase(char ch) {
    if(ch >= 'a' && ch <= 'z') {
        return ch - ('a' - 'A');
    } else {
        return ch;
    }
}

char hexToChar(uint64_t hex) {
    switch(hex) {
    case HEX_FOR_LETTER_A:
        return 'A';
    case HEX_FOR_LETTER_C:
        return 'C';
    case HEX_FOR_LETTER_G:
        return 'G';
    case HEX_FOR_LETTER_T:
        return 'T';
    default:
        return '*';
    }
}

uint64_t charToHex(char ch) {
    ch = UpperCase(ch);
    switch(ch) {
    case 'A':
        return HEX_FOR_LETTER_A;
    case 'C':
        return HEX_FOR_LETTER_C;
    case 'G':
        return HEX_FOR_LETTER_G;
    case 'T':
        return HEX_FOR_LETTER_T;
    default:
        return 0x0;
    }
}

uint64_t getInverseBaseHex(uint64_t base) {
    uint64_t fullHex = (1 << (64 / CHAR_NUM_PER_HEX)) - 1;
    return (~base) & fullHex;
}

uint64_t min_uint64_t(uint64_t value1, uint64_t value2) {
    if(value1 < value2) {
        return value1;
    } else {
        return value2;
    }
}

uint64_t max_uint64_t(uint64_t value1, uint64_t value2) {
    if(value1 > value2) {
        return value1;
    } else {
        return value2;
    }
}


/*
 * Static functions. (file-localized functions)
 */




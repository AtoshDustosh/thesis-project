#include "AuxiliaryDataType.h"

#include <string.h>


static void _StringBufferTest();
static void _HexCodedStringBufferTest();
static void _compareHexCodedStringBufferTest();

void _AuxiliaryDataTypeTestSet() {
    _StringBufferTest();
    _HexCodedStringBufferTest();
    _compareHexCodedStringBufferTest();
}



/**
 * Test the data type StringBuffer.
 */
static void _StringBufferTest() {
    printf("\n**************** _StringBufferTest ****************\n");
    StringBuffer* strBuf1 = NULL;
    StringBuffer* strBuf2 = NULL;

    char* str1 = "1234567";
    char* str2 = "1234";
    strBuf1 = constructStringBuffer(str1, strlen(str1));
    strBuf2 = constructStringBuffer(str2, strlen(str2));

    printStringBuffer(strBuf1);
    printStringBuffer(strBuf2);

    printf("... do string copy from string buffer 1 to string buffer 2\n");
//    strcpy(strBuf2->buffer, strBuf1->buffer); // this will result in an error
    strBuf2->buffer = strBuf1->buffer;

    printStringBuffer(strBuf1);
    printStringBuffer(strBuf2);

    clearStringBuffer(strBuf1);
    clearStringBuffer(strBuf2);
}

/**
 * Test the data type HexCodedStringBuffer.
 */
static void _HexCodedStringBufferTest() {
    printf("\n**************** _HexCodedStringBufferTest ****************\n");
    HexCodedStringBuffer* hexCodedStrBuf1 = NULL;
    HexCodedStringBuffer* hexCodedStrBuf2 = NULL;

    uint64_t arrayForTest1[] = {0x1321adf1353, 0x1321adf1353};
    uint64_t arrayLength1 = 2;
    uint64_t strLength1 = 34;
    uint64_t arrayForTest2[] = {0x0123456789abcdef};
    uint64_t arrayLength2 = 1;
    uint64_t strLength2 = 30;

    hexCodedStrBuf1 = constructHexCodedStringBuffer(arrayForTest1, arrayLength1, strLength1);
    hexCodedStrBuf2 = constructHexCodedStringBuffer(arrayForTest2, arrayLength2, strLength2);

    printHexCodedStringBuffer(hexCodedStrBuf1);
    printHexCodedStringBuffer(hexCodedStrBuf2);

    clearHexCodedStringBuffer(hexCodedStrBuf1);
    clearHexCodedStringBuffer(hexCodedStrBuf2);
}

/**
 * Test function compareHexCodedStringBuffer.
 */
static void _compareHexCodedStringBufferTest() {
    printf("\n**************** _compareHexCodedStringBufferTest ****************\n");
    HexCodedStringBuffer* hexCodedStrBuf1 = NULL;
    HexCodedStringBuffer* hexCodedStrBuf2 = NULL;

    uint64_t arrayForTest1[] = {0x1321adf1353, 0x1321adf1353};
    uint64_t arrayLength1 = 2;
    uint64_t strLength1 = 34;
    uint64_t arrayForTest2[] = {0x0123456789abcdef};
    uint64_t arrayLength2 = 1;
    uint64_t strLength2 = 30;

    printf("... construct hex-coded string buffer 1\n");
    hexCodedStrBuf1 = constructHexCodedStringBuffer(arrayForTest1, arrayLength1, strLength1);
    printHexCodedStringBuffer(hexCodedStrBuf1);

    printf("... construct hex-coded string buffer 2\n");
    hexCodedStrBuf2 = constructHexCodedStringBuffer(arrayForTest2, arrayLength2, strLength2);
    printHexCodedStringBuffer(hexCodedStrBuf2);

    printf("\n");

    printf("compare hex-buffer-1 and hex-buffer-2: %"PRIu64"\n",
           compareHexCodedStringBuffer(hexCodedStrBuf1, hexCodedStrBuf2));
    printf("compare hex-buffer-1 and hex-buffer-1: %"PRIu64"\n",
           compareHexCodedStringBuffer(hexCodedStrBuf1, hexCodedStrBuf1));
    printf("compare hex-buffer-2 and hex-buffer-2: %"PRIu64"\n",
           compareHexCodedStringBuffer(hexCodedStrBuf2, hexCodedStrBuf2));

    clearHexCodedStringBuffer(hexCodedStrBuf1);
    clearHexCodedStringBuffer(hexCodedStrBuf2);
}



/*
 * Working functions.
 */


HexCodedStringBuffer* constructHexCodedStringBuffer(uint64_t* hexArray,
              uint64_t arrayLength, uint64_t strLength) {
    HexCodedStringBuffer* hexCodedStrBuf =
        (HexCodedStringBuffer*)malloc(sizeof(HexCodedStringBuffer));
    if(hexCodedStrBuf == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    if(hexArray == NULL || arrayLength == 0) {
        hexCodedStrBuf->hexArray = NULL;
    } else {
        hexCodedStrBuf->hexArray = (uint64_t*)malloc(sizeof(uint64_t) * arrayLength);
        if(hexCodedStrBuf->hexArray == NULL) {
            printf("ERROR: System memory not enough. \n");
            exit(EXIT_FAILURE);
        }
        for(uint64_t i = 0; i < arrayLength; i++) {
            hexCodedStrBuf->hexArray[i] = hexArray[i];
        }
    }

    hexCodedStrBuf->arrayLength = arrayLength;
    hexCodedStrBuf->strLength = strLength;
    return hexCodedStrBuf;
}

StringBuffer* constructStringBuffer(char* buffer, uint64_t length) {
    StringBuffer* strBuf = (StringBuffer*)malloc(sizeof(StringBuffer));
    if(strBuf == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    strBuf->buffer = (char*)malloc(sizeof(char) * (length + 1));
    if(strBuf->buffer == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    for(uint64_t i = 0; i < length; i++) {
        strBuf->buffer[i] = buffer[i];
    }

    strBuf->buffer[length] = '\0';
    strBuf->length = length;
    return strBuf;
}


uint64_t compareHexCodedStringBuffer(HexCodedStringBuffer* hexCodedStrBuf1,
                                     HexCodedStringBuffer* hexCodedStrBuf2) {
    if(hexCodedStrBuf1->strLength != hexCodedStrBuf2->strLength ||
            hexCodedStrBuf1->arrayLength != hexCodedStrBuf2->arrayLength) {
//        printf("string length or array length not equal. \n");
        return HEX_CODED_STRINGBUFFER_DIFFERNET;
    }
    uint64_t arrayLength = hexCodedStrBuf1->arrayLength;
    for(uint64_t i = 0; i < arrayLength; i++) {
        uint64_t hexInt1 = hexCodedStrBuf1->hexArray[i];
        uint64_t hexInt2 = hexCodedStrBuf2->hexArray[i];

        uint64_t cmpValue = hexInt1 ^ hexInt2;
//        printf("cmpValue: %#"PRIx64"\n", cmpValue);
        /** < \note calculation priority: "!=" > "^". */
        if(cmpValue != 0) {
//            printf("unequal: %#"PRIx64", %#"PRIx64" -> %#"PRIx64"\n",
//                   hexInt1, hexInt2, cmpValue);
            return HEX_CODED_STRINGBUFFER_DIFFERNET;
        }
    }
    return HEX_CODED_STRINGBUFFER_SAME;
}

void printStringBuffer(StringBuffer* strBuf) {
    if(strBuf == NULL) {
        printf("ERROR: null pointer occurred when printing a string buffer. \n");
        exit(EXIT_FAILURE);
    }
    printf("string buffer (0x%p) - (%s, %"PRIu64")\n",
           strBuf, strBuf->buffer, strBuf->length);
}

void printHexCodedStringBuffer(HexCodedStringBuffer* hexCodedStrBuf) {
    if(hexCodedStrBuf == NULL) {
        printf("ERROR: null pointer occurred when printing a hex-coded string buffer. \n");
        exit(EXIT_FAILURE);
    }
    uint64_t i = 0;
    printf("hex-coded string buffer (0x%p) - ({", hexCodedStrBuf);
    uint64_t arrayLength = hexCodedStrBuf->arrayLength;
    for(i = 0; i < arrayLength; i++) {
        printf("%#16"PRIx64"", hexCodedStrBuf->hexArray[i]);
        if(i != hexCodedStrBuf->arrayLength - 1) {
            printf(",");
        }
    }
    printf("}, %"PRIu64", %"PRIu64")\n",
           hexCodedStrBuf->arrayLength, hexCodedStrBuf->strLength);
}

void clearStringBuffer(StringBuffer* strBuf) {
    if(strBuf == NULL) {
        printf("ERROR: null pointer occurs when clearing a string buffer. \n");
        exit(EXIT_FAILURE);
    }
    free(strBuf->buffer);
    free(strBuf);
}

void clearHexCodedStringBuffer(HexCodedStringBuffer* hexCodedStrBuf) {
    if(hexCodedStrBuf == NULL) {
        printf("ERROR: null pointer occurs when clearing a hex-coded string buffer. \n");
        exit(EXIT_FAILURE);
    }
    free(hexCodedStrBuf->hexArray);
    free(hexCodedStrBuf);
}





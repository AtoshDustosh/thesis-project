#include "FileOperation.h"

#include <string.h>

#include "ArrayOperation.h"
#include "AuxiliaryFunction.h"
#include "MyArgs.h"


static void loadOneReadFromFile_projectBufferToRead(char* buffer, uint64_t fastqLine, Read* read);


uint64_t fnaDataSize(char* filePath) {
    printf("\nCalculating size of file %s ... \n", filePath);
    FILE* fp = fopen(filePath, "r");
    uint64_t dataZone = 0;   // if the file pointer is now in the data zone
    uint64_t dataLength = 0;

    if(fp != NULL) {
        char ch = fgetc(fp);
        while(ch != EOF) {
            if(ch == '\n' && !dataZone) {
                // according to format of *.fna file
                // data part is after the first line
                dataZone = 1;
            }
            if(dataZone && ch != '\n') {
                // if not '\n' and is already in the data zone
                dataLength++;
            }
            ch = fgetc(fp);
        }
    } else {
        printf("failed to open file %s", filePath);
    }
    free(fp);

    printf("data size in file %s is %"PRIu64" bp\n", filePath, dataLength);
    return dataLength;
}

void loadFnaData(char* filePath, uint64_t dataLength, uint64_t* hexCodedDNA, char* fnaFileHeader) {
    printf("\nLoading data from file %s ...\n", filePath);
    const uint64_t charNumPerHex = CHAR_NUM_PER_HEX;

    FILE* fp = fopen(filePath, "r");
    if(fp == NULL) {
        printf("failed to open file %s\n", filePath);
        exit(EXIT_FAILURE);
    }

    char ch = fgetc(fp);
    StringBuffer* strBuf = NULL;
    HexCodedStringBuffer* hexCodedStrBuf = NULL;
    uint64_t i = 0;

    char* buffer = (char*)malloc(sizeof(char) * (charNumPerHex + 1));
    if(buffer == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }

    uint64_t hexInt = 0;
    uint64_t strLength = 0;
    uint64_t dataZone = 0;   // if the file pointer is now in the data zone
    uint64_t fileHeaderPointer = 0;
    while(ch != EOF && i < dataLength) {
        // according to format of *.fna file, data part is after the first line
        if(dataZone == 0) {
            if(ch == '\n') {
                dataZone = 1;
                *(fnaFileHeader + fileHeaderPointer) = '\0';
                ch = fgetc(fp);
                continue;
            } else {
                *(fnaFileHeader + fileHeaderPointer) = ch;
                fileHeaderPointer++;
                ch = fgetc(fp);
                continue;
            }
        } else {
            if(ch == '\n') {
                ch = fgetc(fp);
                continue;
            }
        }

        buffer[strLength++] = ch;
        if(strLength % charNumPerHex == 0) {
            /**
             * \alert extract a function if you have the time and fix this pile of shit.
             */
            buffer[strLength] = '\0';
            strBuf = constructStringBuffer(buffer, charNumPerHex);
            hexCodedStrBuf = transStringBufferToHexCodedStringBuffer(strBuf, charNumPerHex);
            hexInt = hexCodedStrBuf->hexArray[0];
            hexCodedDNA[i++] = hexInt;
            clearStringBuffer(strBuf);
            clearHexCodedStringBuffer(hexCodedStrBuf);
//            printf("%s\t", buffer);
//            printf("0x%16"PRIx64"\n", hexInt);
//            if(i % 4 == 0) {
//                printf("\n");
//            }
            buffer = (char*)malloc(sizeof(char) * (charNumPerHex + 1));
            if(buffer == NULL) {
                printf("ERROR: System memory not enough. \n");
                exit(EXIT_FAILURE);
            }
            strLength = 0;
        }
        ch = fgetc(fp);
    }
    buffer[strLength] = '\0';
    strBuf = constructStringBuffer(buffer, charNumPerHex);
    hexCodedStrBuf = transStringBufferToHexCodedStringBuffer(strBuf, charNumPerHex);
    hexInt = hexCodedStrBuf->hexArray[0];
    hexCodedDNA[i++] = hexInt;
    clearStringBuffer(strBuf);
    clearHexCodedStringBuffer(hexCodedStrBuf);
//    printf("%s\t", buffer);
//    printf("0x%16"PRIx64"\n", hexInt);

    free(fp);
}

uint64_t loadOneReadFromFile(char* filePath, FILE** fpointer, Read* read) {
    printf("\n");
    if(*fpointer == NULL) {
        printf("Open file %s\n", filePath);
        *fpointer = fopen(filePath, "r");
    }

    if(*fpointer == NULL) {
        printf("failed to open file %s", filePath);
        exit(EXIT_FAILURE);
    }

    uint64_t bufPointer = 0;
    char buffer[BUFSIZ];
    clearCharArray(buffer, BUFSIZ);

    uint64_t fastqLine = 0;

    printf("Load a read form file %s ... \n", filePath);
    char ch = ' ';
    while(ch != EOF && fastqLine < 4) {
        ch = fgetc(*fpointer);
        if(ch == '\n') {
            buffer[bufPointer] = '\0';
            bufPointer = 0;
            loadOneReadFromFile_projectBufferToRead(buffer, fastqLine, read);
            clearCharArray(buffer, BUFSIZ);
            fastqLine++;
            continue;
        } else if (ch == EOF) {
            break;
        }
        switch(fastqLine) {
        case 0: // header line
            buffer[bufPointer++] = ch;
            break;
        case 1: // sequence line
            buffer[bufPointer++] = UpperCase(ch);
            break;
        case 2: // supplementary line
            break;
        case 3: // quality line
            buffer[bufPointer++] = ch;
            break;
        default:
            printf("ERROR: error occurs when loading a read. \n");
            exit(EXIT_FAILURE);
        }
    }
    if(ch == EOF){
        printf("... reaches end of the fastq file. \n");
        return 0;
    } else {
        return 1;
    }

}




/*
 * Static functions. (file-localized functions)
 */

/**
 * Project the information in buffer to struct type Read.
 *
 * @param buffer buffer that stores information loaded from ?.fastq file
 * @param fastqLine line index of a read's information in ?.fastq file
 * @param read a read
 */
static void loadOneReadFromFile_projectBufferToRead(char* buffer, uint64_t fastqLine, Read* read) {
//    printf("fastq line: %"PRIu64"\n", fastqLine);

    switch(fastqLine) {
    case 0:
        strcpy(read->QNAME, buffer);
        break;
    case 1:
        strcpy(read->SEQ, buffer);
        break;
    case 2:
        break;
    case 3:
        strcpy(read->QUAL, buffer);
        break;
    default:
        printf("Function loadOneReadFromFile_projectBufferToRead args error. \n");
        exit(EXIT_FAILURE);
    }
}
























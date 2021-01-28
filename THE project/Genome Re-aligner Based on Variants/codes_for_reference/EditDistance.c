#include "EditDistance.h"

#include <string.h>
#include <strings.h>

#include "AuxiliaryFunction.h"
#include "Queue.h"

#define DIAGONALLYEXTENDED 1

#define INITEDVALUE 9999


/*
 * Fill the edit-distance matrix.
 */
static void fillEditDistanceMatrix(StringBuffer* strRow, StringBuffer* strColumn,
                                   uint64_t* EDmatrix, uint64_t* diagonallyExtendedMatrix,
                                   const uint64_t EDmax);
// auxiliary functions - extracted from detailed steps
static void extendEDMatrixDiagonally(StringBuffer* strRow, StringBuffer* strColumn,
                                     uint64_t startRow, uint64_t startColumn, uint64_t* endRow,
                                     uint64_t* endColumn, uint64_t* EDmatrix,
                                     uint64_t* diagonallyExtendedMatrix, uint64_t EDmax);
static void extendEDMatrixHorizontally(StringBuffer* strRow, StringBuffer* strColumn,
                                       uint64_t startRow, uint64_t startColumn, uint64_t* EDmatrix,
                                       uint64_t* diagonallyExtendedMatrix, uint64_t EDmax);
static void extendEDMatrixVertically(StringBuffer* strRow, StringBuffer* strColumn,
                                     uint64_t startRow, uint64_t startColumn, uint64_t* EDmatrix,
                                     uint64_t* diagonallyExtendedMatrix, uint64_t EDmax);
static void splitAfterExtension(StringBuffer* strRow, StringBuffer* strColumn,
                                Queue* startRowQueue, Queue*startColumnQueue,
                                uint64_t endRow, uint64_t endColumn, uint64_t* EDmatrix,
                                uint64_t* diagonallyExtendedMatrix, uint64_t EDmax);
/*
 * Process the edit-distance matrix and get the smallest edit-distance and its best string interval
 * in strColumn (the ref string) that contains a pattern of strRow (pattern string).
 */
static uint64_t processEDMatrix(StringBuffer* strRow, StringBuffer* strColumn, uint64_t* EDmatrix,
                                char* CIGARbuffer, uint64_t maxBufLen, const uint64_t EDmax);
static void parseCIGAR(char* CIGARbuffer);

/*
 * Declarations of tests.
 */
static void _calculateEditDistanceTest();


void _EditDistanceTestSet() {
    _calculateEditDistanceTest();
}


/*
 * Tests for working functions.
 */

/**
 * Test function calculateEditDistance.
 */
static void _calculateEditDistanceTest() {
    printf("\n**************** _calculateEditDistanceTest ****************\n");
    StringBuffer* strBufRow = NULL;
    StringBuffer* strBufColumn = NULL;
    uint64_t maxBufLen = BUFSIZ;
    char* CIGAR = (char*)malloc(sizeof(char) * maxBufLen);
    if(CIGAR == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }

    uint64_t EDmax = 4;
    uint64_t bestED = 0;

    char* string1 = NULL;
    char* string2 = NULL;

    string1 = "agtcgccgctgctgc";
    string2 = "agcgcttgctgc";
//    string1 = "agtcgccgctgctgc";  // out of EDmax
//    string2 = "agggggggtgc";
    strBufColumn = constructStringBuffer(string1, (uint64_t)strlen(string1));
    strBufRow = constructStringBuffer(string2, (uint64_t)strlen(string2));
    bestED = calculateEditDistance(strBufRow, strBufColumn, EDmax, CIGAR, maxBufLen);
    printf("Best edit-distance: %"PRIu64"\n", bestED);
    printf("CIGAR: %s\n", CIGAR);
    clearStringBuffer(strBufColumn);
    clearStringBuffer(strBufRow);

    printf("******************************************************************\n");
    string1 = "ccagtcgctgcgctt";
    string2 = "agcgcttgcgc";
    strBufColumn = constructStringBuffer(string1, (uint64_t)strlen(string1));
    strBufRow = constructStringBuffer(string2, (uint64_t)strlen(string2));
    bestED = calculateEditDistance(strBufRow, strBufColumn, EDmax, CIGAR, maxBufLen);
    printf("Best edit-distance: %"PRIu64"\n", bestED);
    printf("CIGAR: %s\n", CIGAR);
    clearStringBuffer(strBufColumn);
    clearStringBuffer(strBufRow);

    printf("******************************************************************\n");
    string1 = "ccagtcgctgcgcttacac";
    string2 = "ccagtcgctgcgctt";
    strBufColumn = constructStringBuffer(string1, (uint64_t)strlen(string1));
    strBufRow = constructStringBuffer(string2, (uint64_t)strlen(string2));
    bestED = calculateEditDistance(strBufRow, strBufColumn, EDmax, CIGAR, maxBufLen);
    printf("Best edit-distance: %"PRIu64"\n", bestED);
    printf("CIGAR: %s\n", CIGAR);
    clearStringBuffer(strBufColumn);
    clearStringBuffer(strBufRow);

    printf("******************************************************************\n");
    string1 = "ttttccagtcgctgcgctt";
    string2 = "ccagtcgctgcgctt";
    strBufColumn = constructStringBuffer(string1, (uint64_t)strlen(string1));
    strBufRow = constructStringBuffer(string2, (uint64_t)strlen(string2));
    bestED = calculateEditDistance(strBufRow, strBufColumn, EDmax, CIGAR, maxBufLen);
    printf("Best edit-distance: %"PRIu64"\n", bestED);
    printf("CIGAR: %s\n", CIGAR);
    clearStringBuffer(strBufColumn);
    clearStringBuffer(strBufRow);

    printf("******************************************************************\n");
    string1 = "tttttccagtcgctgcgctt";
    string2 = "ccagtcgctgcgctt";
    strBufColumn = constructStringBuffer(string1, (uint64_t)strlen(string1));
    strBufRow = constructStringBuffer(string2, (uint64_t)strlen(string2));
    bestED = calculateEditDistance(strBufRow, strBufColumn, EDmax, CIGAR, maxBufLen);
    printf("Best edit-distance: %"PRIu64"\n", bestED);
    printf("CIGAR: %s\n", CIGAR);
    clearStringBuffer(strBufColumn);
    clearStringBuffer(strBufRow);


    free(CIGAR);
}

/*
 * Working functions.
 */
uint64_t calculateEditDistance(StringBuffer* patternStrBuf, StringBuffer* refStrBuf, uint64_t EDmax,
                               char* CIGARbuffer, uint64_t maxBufLen) {
    /*
     * Actually, considering the size of EDmax, we don't need to use uint64_t, and uint8_t
     * is enough for EDmax and scoreMatrix.
     * But for better correctness, we leave optimizations for future work.
     */
    const uint64_t rowNum = patternStrBuf->length + 1;
    const uint64_t columnNum = refStrBuf->length + 1;

    uint64_t bestED = 0;

    /*
     * Initializations.
     */
    char* strRow = (char*)malloc(sizeof(char) * (rowNum + 1));
    char* strColumn = (char*)malloc(sizeof(char) * (columnNum + 1));
    if(strRow == NULL || strColumn == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    uint64_t EDmatrix[rowNum][columnNum];
    uint64_t diagonallyExtendedMatrix[rowNum][columnNum];

    for(uint64_t i = 0; i < rowNum; i++) {      // initialize edit-distance matrix
        for(uint64_t j = 0; j < columnNum; j++) {
            EDmatrix[i][j] = INITEDVALUE;
            diagonallyExtendedMatrix[i][j] = 0;
        }
    }
    for(uint64_t i = 0; i < rowNum;
            i++) {      // initialize first row (bounder) of edit-distance matrix
        EDmatrix[i][0] = i;
    }
    for(uint64_t i = 0; i < columnNum;
            i++) {      // initialize first column (bounder) of edit-distance matrix
        EDmatrix[0][i] = i;
    }
    for(uint64_t i = 0; i < rowNum; i++) {      // copy patternStrBuf
        strRow[i + 1] = patternStrBuf->buffer[i];
    }
    strRow[0] = ' ';
    strRow[rowNum] = '\0';
    for(uint64_t i = 0; i < columnNum; i++) {   // copy refStrBuf
        strColumn[i + 1] = refStrBuf->buffer[i];
    }
    strColumn[0] = ' ';
    strColumn[columnNum] = '\0';

    StringBuffer* strBufRow = constructStringBuffer(strRow, rowNum);
    StringBuffer* strBufColumn = constructStringBuffer(strColumn, columnNum);

//    printf("reference: ");
//    printStringBuffer(strBufColumn);
//    printf("pattern:   ");
//    printStringBuffer(strBufRow);
//    printf("ED matrix - rowNum: %"PRIu64", columnNum: %"PRIu64"\n",
//           strBufRow->length, strBufColumn->length);

    /*
     * Calculate the edit-distance matrix.
     */
//    printf("initialized ED matrix:\n");
//    for(uint64_t i = 0; i < rowNum; i++) {
//        for(uint64_t j = 0; j < columnNum; j++) {
//            printf("%"PRIu64"\t", EDmatrix[i][j]);
//        }
//        printf("\n\n");
//    }

//    printf("rowNum: %"PRIu64", columnNum: %"PRIu64"\n", rowNum, columnNum);
    fillEditDistanceMatrix(strBufRow, strBufColumn, (uint64_t*)EDmatrix,
                           (uint64_t*)diagonallyExtendedMatrix, EDmax);


//    printf("calculated ED matrix:\n");
//    for(uint64_t i = 0; i < rowNum; i++) {
//        for(uint64_t j = 0; j < columnNum; j++) {
//            printf("%"PRIu64"\t", EDmatrix[i][j]);
//        }
//        printf("\n\n");
//    }

    /**< \todo handle CIGAR string - generated from EDmatrix */
    bestED = processEDMatrix(strBufRow, strBufColumn, (uint64_t*)EDmatrix, CIGARbuffer, maxBufLen,
                             EDmax);
//    printf(" in calculation of edit-distance - print CIGAR: %s\n", CIGARbuffer);

    clearStringBuffer(strBufRow);
    clearStringBuffer(strBufColumn);
    return bestED;
}


/*
 * Static functions. (file-localized functions)
 */

/**
 * Fill the edit-distance matrix of 2 strings until filling-work cannot be done any more.
 *
 * @param strRow row string of the edit-distance matrix
 * @param strColumn column string of the edit-distance matrix
 * @param EDmatrix edit-distance matrix - for i = 0 : n-1 {EDmatrix[0][i] = EDmatrix[i][0] = i};
 *      all other values are initialized as INITEDVALUE
 * @param diagonallyExtendedMatrix flags marking whether a point has been tried to extend diagonally
 * @param EDmax maximum edit distance that can be allowed
 */
static void fillEditDistanceMatrix(StringBuffer* strRow, StringBuffer* strColumn,
                                   uint64_t* EDmatrix, uint64_t* diagonallyExtendedMatrix,
                                   const uint64_t EDmax) {
    /**
     * \note for "matrix[i][j]", use "matrix[row * columnNum + column]",
     *       or "*(matrix + row * columnNum + column)".
     */

//    const uint64_t rowNum = strRow->length;
//    const uint64_t columnNum = strColumn->length;

    uint64_t startRow = 0;
    uint64_t startColumn = 0;

    uint64_t endRow = 0;
    uint64_t endColumn = 0;

    Queue* startRowQueue = (Queue*)malloc(sizeof(Queue));
    Queue* startColumnQueue = (Queue*)malloc(sizeof(Queue));
    QueueCell* queueCell = (QueueCell*)malloc(sizeof(QueueCell));
    if(queueCell == NULL || startColumnQueue == NULL || startRowQueue == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }

    /**< \note initializing enqueue operations */
    startRowQueue = initQueue();
    startColumnQueue = initQueue();
    enQueue(startRowQueue, newQueueCell(startRow));
    enQueue(startColumnQueue, newQueueCell(startColumn));

    uint64_t loopCount = 0;
    uint64_t row = 0;
    uint64_t column = 0;
    while(startRowQueue->length != 0 && startColumnQueue->length != 0) {
        loopCount++;
        /**< \note dequeue operation */
        deQueue(startRowQueue, queueCell);
        startRow = queueCell->data;
        deQueue(startColumnQueue, queueCell);
        startColumn = queueCell->data;

//        printf(">>>>> dequeue and start from (%"PRIu64", %"PRIu64")\n", startRow, startColumn);
//        printStringBuffer(strRow);
//        printStringBuffer(strColumn);
//        printf("ED-matrix:\n");
//        for(uint64_t m = 0; m < rowNum; m++) {
//            for(uint64_t n = 0; n < columnNum; n++) {
//                printf("%"PRIu64"\t", *(EDmatrix + m * columnNum + n));
//            }
//            printf("\n\n");
//        }
//        printf("Diagonally-extended-matrix:\n");
//        for(uint64_t m = 0; m < rowNum; m++) {
//            for(uint64_t n = 0; n < columnNum; n++) {
//                printf("%"PRIu64"\t", *(diagonallyExtendedMatrix + m * columnNum + n));
//            }
//            printf("\n\n");
//        }

        endRow = startRow;
        endColumn = startColumn;

        /**< \note extend diagonally until encountering a mismatch, reaching ends of the matrix,
            or encountering a point that has already been extended diagonally*/
        extendEDMatrixDiagonally(strRow, strColumn, startRow, startColumn, &endRow, &endColumn,
                                 EDmatrix, diagonallyExtendedMatrix, EDmax);
//        printf("extend diagonally:(%"PRIu64",%"PRIu64") to (%"PRIu64",%"PRIu64")\n",
//               startRow, startColumn, endRow, endColumn);
//        if(endRow == rowNum || endColumn == columnNum) {
//            printf("ERROR: endRow and endColumn out of range. \n");
//            exit(EXIT_FAILURE);
//        }

        /**< \note extend horizontally and vertically from all points that have just been extended
        diagonally */
        row = startRow;
        column = startColumn;
        for(; row <= endRow && column <= endColumn; row++, column++) {
            /**< \note extend horizontally within the range of  "ED <= EDmax" */
            extendEDMatrixHorizontally(strRow, strColumn, row, column + 1, EDmatrix,
                                       diagonallyExtendedMatrix, EDmax);
            /**< \note extend vertically within the range of "ED <= EDmax*/
            extendEDMatrixVertically(strRow, strColumn, row + 1, column, EDmatrix,
                                     diagonallyExtendedMatrix, EDmax);
        }

        /**< \note process splitting operations */
        splitAfterExtension(strRow, strColumn, startRowQueue, startColumnQueue, endRow, endColumn,
                            EDmatrix, diagonallyExtendedMatrix, EDmax);
    }

//    printf(">>>>>>>>>> total loop count: %"PRIu64" >>>>> ", loopCount);
//    printf("ED-matrix:\n");
//    for(uint64_t m = 0; m < rowNum; m++) {
//        for(uint64_t n = 0; n < columnNum; n++) {
//            printf("%"PRIu64"\t", *(EDmatrix + m * columnNum + n));
//        }
//        printf("\n\n");
//    }
//    printf("Diagonally-extended-matrix:\n");
//    for(uint64_t m = 0; m < rowNum; m++) {
//        for(uint64_t n = 0; n < columnNum; n++) {
//            printf("%"PRIu64"\t", *(diagonallyExtendedMatrix + m * columnNum + n));
//        }
//        printf("\n\n");
//    }


    clearQueue(startRowQueue);
    clearQueue(startColumnQueue);
    free(queueCell);
}

static void extendEDMatrixDiagonally(StringBuffer* strRow, StringBuffer* strColumn,
                                     uint64_t startRow, uint64_t startColumn, uint64_t* endRow,
                                     uint64_t* endColumn, uint64_t* EDmatrix,
                                     uint64_t* diagonallyExtendedMatrix, uint64_t EDmax) {
    uint64_t rowNum = strRow->length;
    uint64_t columnNum = strColumn->length;

    uint64_t row = startRow;
    uint64_t column = startColumn;

    *endRow = startRow;
    *endColumn = startColumn;
    for(; row < rowNum && column < columnNum; row++, column++) {
        if(row == 0 || column == 0) {
            continue;
        }
//        if(*endRow == rowNum - 1 || *endColumn == columnNum - 1) {
//            break;
//        }
        char charOfThisRow = strRow->buffer[row];
        char charOfThisColumn = strColumn->buffer[column];
        uint64_t EDofLastPoint = *(EDmatrix + (row - 1) * columnNum + (column - 1));
        uint64_t EDofLeftPoint = *(EDmatrix + row * columnNum + (column - 1));
        uint64_t EDofUpperPoint =  *(EDmatrix + (row - 1) * columnNum + column);

        uint64_t newEDvalue = 0;
//        printf("(%"PRIu64", %"PRIu64") (%c, %c)\n", row, column, charOfThisRow, charOfThisColumn);
        if(charOfThisRow == charOfThisColumn) {
            newEDvalue = EDofLastPoint;
            if(newEDvalue > EDmax) {
                break;
            }
            *(EDmatrix + row * columnNum + column) = newEDvalue;
            *(diagonallyExtendedMatrix + row * columnNum + column) = DIAGONALLYEXTENDED;
            *endRow = *endRow + 1;
            *endColumn = *endColumn + 1;
        } else if(charOfThisRow != charOfThisColumn
                  && DIAGONALLYEXTENDED != *(diagonallyExtendedMatrix + row * columnNum + column)
                 ) {
            newEDvalue = min_uint64_t(min_uint64_t(EDofLastPoint, EDofLeftPoint),
                                      min_uint64_t(EDofLastPoint, EDofUpperPoint)) + 1;
            if(newEDvalue > EDmax) {
                break;
            }
            *(EDmatrix + row * columnNum + column) = newEDvalue;
            /**< \alert don't mark this point as diagonally-extended */
            break;
        } else {
            break;
        }
    }
}

static void extendEDMatrixHorizontally(StringBuffer* strRow, StringBuffer* strColumn,
                                       uint64_t startRow, uint64_t startColumn, uint64_t* EDmatrix,
                                       uint64_t* diagonallyExtendedMatrix, uint64_t EDmax) {
    uint64_t columnNum = strColumn->length;

    for(uint64_t column = startColumn; column < columnNum; column++) {
        char charOfThisRow = strRow->buffer[startRow];
        char charOfThisColumn = strColumn->buffer[column];

        uint64_t EDofLastPoint = *(EDmatrix + (startRow - 1) * columnNum + (column - 1));
        uint64_t EDofLeftPoint = *(EDmatrix + startRow * columnNum + (column - 1));
        uint64_t EDofUpperPoint = *(EDmatrix + (startRow - 1) * columnNum + column);

        if(*(EDmatrix + startRow * columnNum + column) != INITEDVALUE) {
            break;  /**< \note avoid adding the same matrix unit repeatedly */
        }

        uint64_t newEDvalue = 0;
        if(charOfThisRow == charOfThisColumn) {
            newEDvalue = EDofLastPoint;
            if(newEDvalue > EDmax) {
                break;
            }
            *(EDmatrix + startRow * columnNum + column) = newEDvalue;
        } else {
            newEDvalue = min_uint64_t(min_uint64_t(EDofLastPoint, EDofLeftPoint),
                                      min_uint64_t(EDofLastPoint, EDofUpperPoint)) + 1;
            if(newEDvalue > EDmax) {
                break;
            }
            *(EDmatrix + startRow * columnNum + column) = newEDvalue;
        }
    }
}

static void extendEDMatrixVertically(StringBuffer* strRow, StringBuffer* strColumn,
                                     uint64_t startRow, uint64_t startColumn, uint64_t* EDmatrix,
                                     uint64_t* diagonallyExtendedMatrix, uint64_t EDmax) {
    uint64_t rowNum = strRow->length;
    uint64_t columnNum = strColumn->length;

    for(uint64_t row = startRow; row < rowNum; row++) {
        char charOfThisRow = strRow->buffer[row];
        char charOfThisColumn = strColumn->buffer[startColumn];

        uint64_t EDofLastPoint = *(EDmatrix + (row - 1) * columnNum + (startColumn - 1));
        uint64_t EDofLeftPoint = *(EDmatrix + row * columnNum + (startColumn - 1));
        uint64_t EDofUpperPoint = *(EDmatrix + (row - 1) * columnNum + startColumn);

        if(*(EDmatrix + row * columnNum + startColumn) != INITEDVALUE) {
            break;  /**< \note avoid adding the same matrix unit repeatedly */
        }

        uint64_t newEDvalue = 0;
        if(charOfThisRow == charOfThisColumn) {
            newEDvalue = EDofLastPoint;
            if(newEDvalue > EDmax) {
                break;
            }
            *(EDmatrix + row * columnNum + startColumn) = newEDvalue;
        } else {
            newEDvalue = min_uint64_t(min_uint64_t(EDofLastPoint, EDofLeftPoint),
                                      min_uint64_t(EDofLastPoint, EDofUpperPoint)) + 1;
            if(newEDvalue > EDmax) {
                break;
            }
            *(EDmatrix + row * columnNum + startColumn) = newEDvalue;
        }
    }
}

static void splitAfterExtension(StringBuffer* strRow, StringBuffer* strColumn,
                                Queue* startRowQueue, Queue*startColumnQueue,
                                uint64_t endRow, uint64_t endColumn, uint64_t* EDmatrix,
                                uint64_t* diagonallyExtendedMatrix, uint64_t EDmax) {
    /* conditions:
    1. endRow + 1 == rowNum && endColumn + 1 == columnNum - finishes - end the process
    2. endRow + 1 == rowNum && endColumn + 1 < columnNum - row ends
    3. endRow + 1 < rowNum && endColumn + 1 == columnNum - column ends
    4. endRow + 1 < rowNum && endColumn + 1 < columnNum
        4.1 EDof(endRow, endColumn) != EDof(endRow + 1, endColumn + 1) &&
            EDof(endRow + 1, endColumn + 1) <= EDmax - mismatch
        4.2 EDof(endRow, endColumn) == EDof(endRow + 1, endColumn + 1) &&
            EDof(endRow + 1, endColumn + 1) <= EDmax - match
        4.3 EDof(endRow + 1, endColumn + 1) > EDmax - exceed range of EDmax
    5. ... (NULL)
    */
    /**< \note the order of the following cases are based on
        priority: min ED point > next point > right point > lower point */
    uint64_t rowNum = strRow->length;
    uint64_t columnNum = strColumn->length;

    if(endRow + 1 == rowNum && endColumn + 1 == columnNum) { // procedure finishes
//        printf("reaches end of matrix. \n");
        return;
    } else if(endRow + 1 == rowNum && endColumn + 1 < columnNum) {  // reaches bounder of row
        char charOfThisRow = strRow->buffer[endRow];
        char charOfThisColumn = strColumn->buffer[endColumn + 1];

        uint64_t EDofLastPoint = *(EDmatrix + (endRow - 1) * columnNum + endColumn);
        uint64_t EDofLeftPoint = *(EDmatrix + endRow * columnNum + endColumn);
        uint64_t EDofUpperPoint = *(EDmatrix + (endRow - 1) * columnNum + (endColumn + 1));

        uint64_t newEDvalue = 0;
        if(charOfThisRow == charOfThisColumn) {
            newEDvalue = EDofLastPoint;
            if(newEDvalue > EDmax) {
                return;
            }
            *(EDmatrix + endRow * columnNum + (endColumn + 1)) = newEDvalue;
        } else {
            newEDvalue = min_uint64_t(min_uint64_t(EDofLastPoint, EDofLeftPoint),
                                      min_uint64_t(EDofLastPoint, EDofUpperPoint)) + 1;
            if(newEDvalue > EDmax) {
                return;
            }
            *(EDmatrix + endRow * columnNum + (endColumn + 1)) = newEDvalue;
        }
    } else if(endRow + 1 < rowNum && endColumn + 1 == columnNum) {  // reaches bounder of column
        char charOfThisRow = strRow->buffer[endRow + 1];
        char charOfThisColumn = strColumn->buffer[endColumn];

        uint64_t EDofLastPoint = *(EDmatrix + endRow * columnNum + (endColumn - 1));
        uint64_t EDofUpperPoint = *(EDmatrix + endRow * columnNum + endColumn);
        uint64_t EDofLeftPoint = *(EDmatrix + (endRow + 1) * columnNum + (endColumn - 1));

        uint64_t newEDvalue = 0;
        if(charOfThisRow == charOfThisColumn) {
            newEDvalue = EDofLastPoint;
            if(newEDvalue > EDmax) {
                return;
            }
            *(EDmatrix + (endRow + 1) * columnNum + endColumn) = newEDvalue;
        } else {
            newEDvalue = min_uint64_t(min_uint64_t(EDofLastPoint, EDofLeftPoint),
                                      min_uint64_t(EDofLastPoint, EDofUpperPoint)) + 1;
            if(newEDvalue > EDmax) {
                return;
            }
            *(EDmatrix + (endRow + 1) * columnNum + endColumn) = newEDvalue;
        }
    } else if(endRow + 1 < rowNum && endColumn + 1 < columnNum) {    // split into 3 streams
        uint64_t continueExtension = 0;
        uint64_t EDofNextPoint = *(EDmatrix + (endRow + 1) * columnNum + (endColumn + 1));
        uint64_t EDofRightPoint = *(EDmatrix + endRow * columnNum + (endColumn + 1));
        uint64_t EDofLowerPoint = *(EDmatrix + (endRow + 1) * columnNum + endColumn);

        if(EDofNextPoint <= EDmax && DIAGONALLYEXTENDED !=
                *(diagonallyExtendedMatrix + (endRow + 1) * columnNum + (endColumn + 1))) {
            enQueue(startRowQueue, newQueueCell(endRow + 1));
            enQueue(startColumnQueue, newQueueCell(endColumn + 1));
            *(diagonallyExtendedMatrix + (endRow + 1) * columnNum + (endColumn + 1)) =
                DIAGONALLYEXTENDED;
            continueExtension++;
//            printf("... (next) can extend to next point (%"PRIu64",%"PRIu64")\n",
//                   endRow + 1, endColumn + 1);
        }
        if(EDofRightPoint <= EDmax && DIAGONALLYEXTENDED !=
                *(diagonallyExtendedMatrix + endRow * columnNum + (endColumn + 1))) {
            enQueue(startRowQueue, newQueueCell(endRow));
            enQueue(startColumnQueue, newQueueCell(endColumn + 1));
            *(diagonallyExtendedMatrix + endRow * columnNum + (endColumn + 1)) =
                DIAGONALLYEXTENDED;
            continueExtension++;
//            printf("... (right) can extend to right point (%"PRIu64",%"PRIu64")\n",
//                   endRow, endColumn + 1);
        }
        if(EDofLowerPoint <= EDmax && DIAGONALLYEXTENDED !=
                *(diagonallyExtendedMatrix + (endRow + 1) * columnNum + endColumn)) {
            enQueue(startRowQueue, newQueueCell(endRow + 1));
            enQueue(startColumnQueue, newQueueCell(endColumn));
            *(diagonallyExtendedMatrix + (endRow + 1) * columnNum + endColumn) =
                DIAGONALLYEXTENDED;
            continueExtension++;
//            printf("... (down) can extend to lower point (%"PRIu64",%"PRIu64")\n",
//                   endRow + 1, endColumn);
        }

        if(continueExtension == 0) {
//            printf("... cannot continue extension. \n");
        }
    }

}


/**
 * Process the edit-distance matrix and get the smallest edit-distance and its best string interval
 * in strColumn (the ref string) that contains a pattern of strRow (pattern string).
 *
 * @param strRow row string of the edit-distance matrix
 * @param strColumn column string of the edit-distance matrix
 * @param EDmatrix edit-distance matrix - for i = 0 : n-1 {EDmatrix[0][i] = EDmatrix[i][0] = i};
 *      all other values are initialized as INITEDVALUE
 * @param CIGARbuffer buffer for CIGAR string; "*" if cannot find a match within EDmax
 * @param maxBufLen maximum length of buffer for CIGAR string
 * @param EDmax maximum edit distance that can be allowed
 * @return best edit-distance if it's within limit; INITEDVALUE if out of limit
 */
static uint64_t processEDMatrix(StringBuffer* strRow, StringBuffer* strColumn, uint64_t* EDmatrix,
                                char* CIGARbuffer, uint64_t maxBufLen, const uint64_t EDmax) {
    uint64_t rowNum = strRow->length;
    uint64_t columnNum = strColumn->length;

    /**
     * \note find the best edit-distance and its column in ED matrix.
     */
    uint64_t bestED = INITEDVALUE;
    uint64_t bestEDColumn = 0;
    for(uint64_t i = 1; i < columnNum; i++) {
        uint64_t EDvalue = *(EDmatrix + (rowNum - 1) * columnNum + i);
        if(EDvalue < bestED) {
            bestED = EDvalue;
            bestEDColumn = i;
        }
    }

    if(bestED == INITEDVALUE) {
        /** < cannot find a match within EDmax */
        CIGARbuffer[0] = '*';
        CIGARbuffer[1] = '\0';
        return bestED;
    }

    /**
     * Reconstruct CIGAR string from ED matrix
     */
    uint64_t row = rowNum - 1;
    uint64_t column = bestEDColumn;
    uint64_t CIGARbufferPointer = 0;
    while(row >= 0 && column >= 0 && CIGARbufferPointer < maxBufLen) {
        uint64_t EDofLastPoint = INITEDVALUE;
        uint64_t EDofLeftPoint = INITEDVALUE;
        uint64_t EDofUpperPoint = INITEDVALUE;
        if(row > 0 && column > 0) {
            EDofLastPoint = *(EDmatrix + (row - 1) * columnNum + (column - 1));
        }
        if(column > 0) {
            EDofLeftPoint = *(EDmatrix + row * columnNum + (column - 1));
        }
        if(row > 0) {
            EDofUpperPoint = *(EDmatrix + (row - 1) * columnNum + column);
        }

        uint64_t minEDvalue = min_uint64_t(min_uint64_t(EDofLastPoint, EDofLeftPoint),
                                           min_uint64_t(EDofLastPoint, EDofUpperPoint));
//        printf("%"PRIu64" (row, column): (%"PRIu64",%"PRIu64")\n", CIGARbufferPointer, row, column);
//        printf("last:%"PRIu64",upper:%"PRIu64",left:%"PRIu64"\n",
//               EDofLastPoint, EDofUpperPoint, EDofLeftPoint);
        if(minEDvalue == EDofLastPoint) {
            CIGARbuffer[CIGARbufferPointer++] = 'M';
            row--;
            column--;
        } else if(minEDvalue == EDofLeftPoint && minEDvalue != EDofLastPoint
                  && minEDvalue != EDofUpperPoint) {
            CIGARbuffer[CIGARbufferPointer++] = 'D';
            column--;
        } else {
            CIGARbuffer[CIGARbufferPointer++] = 'I';
            row--;
        }
        if(row == 0 && column == 0) {
            break;
        }
    }
    CIGARbuffer[CIGARbufferPointer] = '\0';
    reverseString(CIGARbuffer);
//    printf("CIGAR: %s\n", CIGARbuffer);
    parseCIGAR(CIGARbuffer);
    return bestED;
}

/**
 * Parse CIGARbuffer from continuous characters to numbers and characters.
 *
 * @param / @return CIGARbuffer buffer for CIGAR string
 */
static void parseCIGAR(char* CIGARbuffer) {
    uint64_t strLength = strlen(CIGARbuffer);
    char* CIGAR = (char*)malloc(sizeof(char) * (strLength + 1));
    if(CIGAR == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    uint64_t CIGARpointer = 0;
    uint64_t CIGARbufferPointer = 0;
    char ch = '\0';
    char prevCh = ch;
    uint64_t countSameCh = 0;

    /** < \note copy CIGAR buffer into another string */
    for(uint64_t i = 0; i < strLength; i++) {
        CIGAR[i] = CIGARbuffer[i];
    }
    CIGAR[strLength] = '\0';
//    printf("copied CIGAR: %s\n", CIGAR);

    countSameCh = 0;
    CIGARpointer = 0;
    CIGARbufferPointer = 0;
    countSameCh = 0;
    char numString[BUFSIZ];
    uint64_t numStrLength = 0;
    ch = CIGAR[CIGARpointer++];
    prevCh = ch;
    countSameCh++;
    for(; CIGARpointer < strLength; CIGARpointer++) {
        ch = CIGAR[CIGARpointer];
        if(prevCh == ch) {
            countSameCh++;
        } else {
            sprintf(numString, "%"PRIu64, countSameCh);
            numStrLength = (uint64_t)strlen(numString);
//            printf("numStringLength: %"PRIu64", numString: \"%s\"\n", numStrLength, numString);
            for(uint64_t i = 0; i < numStrLength; i++) {
                CIGARbuffer[CIGARbufferPointer++] = numString[i];
            }
            CIGARbuffer[CIGARbufferPointer++] = prevCh;
            countSameCh = 1;
        }
        prevCh = ch;
    }
    sprintf(numString, "%"PRIu64, countSameCh);
    numStrLength = (uint64_t)strlen(numString);
//    printf("numStringLength: %"PRIu64", numString: \"%s\"\n", numStrLength, numString);
    for(uint64_t i = 0; i < numStrLength; i++) {
        CIGARbuffer[CIGARbufferPointer++] = numString[i];
    }
    CIGARbuffer[CIGARbufferPointer++] = prevCh;
    CIGARbuffer[CIGARbufferPointer] = '\0';

//    printf("parsed CIGAR: %s\n", CIGARbuffer);
    free(CIGAR);
}











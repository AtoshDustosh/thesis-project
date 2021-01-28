#include "ArrayOperation.h"


static void _clearCharArrayTest();


void _ArrayOperationTestSet() {
    _clearCharArrayTest();
}


/*
 * Tests for working functions.
 */

/**
 * Test function clearCharArray.
 */
static void _clearCharArrayTest() {
    printf("\n**************** _clearCharArrayTest ****************\n");
    int i = 0;
    char charArray[4] = {'a', 'c', 'g', 't'};

    printf("charArray: ");
    for(i = 0; i < 4; i++) {
        printf("%c.", charArray[i]);
    }
    printf(".\n");

    printf("... clear charArray\n");
    clearCharArray(charArray, 4);

    printf("charArray: ");
    for(i = 0; i < 4; i++) {
        printf("%c.", charArray[i]);
    }
    printf(".\n");

}








/*
 * Working functions.
 */

void clearCharArray(char charArray[], int arrayLength) {
    int i = 0;
    for(i = 0; i < arrayLength; i++) {
        charArray[i] = '\0';
    }
}

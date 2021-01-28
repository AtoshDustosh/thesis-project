#include "cDataType.h"


static void _uint64_tTest();
static void _uintTypeCastTest();


void _cDataTypeTestSet() {
    _uint64_tTest();
    _uintTypeCastTest();
}


/*
 * Functions for testing.
 */

/**
 * Test data type "uint64_t".
 */
static void _uint64_tTest() {
    printf("\n**************** _uint64_tTest ****************\n");

    uint64_t num1 = 8914584519061689668;
    uint64_t num1_hex = 0x7BB6F7637ADE5544;
    uint64_t num2 = 8914584519061689668 >> 1;
    uint64_t num2_hex = 0x7BB6F7637ADE5544 >> 1;
    printf("uint64_t num1: %"PRIu64"\n", num1);
    printf("uint64_t num1_hex: %#"PRIx64"\n", num1_hex);
    printf("uint64_t num2: %"PRIu64"\n", num2);
    printf("uint64_t num2_hex: %#"PRIx64"\n", num2_hex);

    printf("num1 ^ num2: %#"PRIx64"\n", num1_hex ^ num2_hex);
    printf("num1 & num2: %#"PRIx64"\n", num1_hex & num2_hex);

    // bit-shifting will not affect the symbol of uint64_t
    printf("(num1 << 1) >> 1: %#"PRIx64"\n", (num1 << 1) >> 1);
}

/**
 * Test the effect of casting operation n uint_? data types.
 */
static void _uintTypeCastTest(){
    printf("\n**************** _uint64_tTest ****************\n");
    uint8_t uint_8bit = -1;
    uint64_t uint_64bit = 0;

    printf("uint8_t: %"PRIu8"\n", uint_8bit);
    printf("uint8_t: %#"PRIx8"\n", uint_8bit);
    printf("uint64_t: %"PRIu64"\n", uint_64bit);
    printf("uint64_t: %#"PRIx64"\n", uint_64bit);

    printf("... cast uint8_t to uint64_t\n");
    uint_64bit = (uint64_t)uint_8bit;


    printf("uint8_t: %"PRIu8"\n", uint_8bit);
    printf("uint8_t: %#"PRIx8"\n", uint_8bit);
    printf("uint64_t: %"PRIu64"\n", uint_64bit);
    printf("uint64_t: %#"PRIx64"\n", uint_64bit);


}


















/*
 * Working functions.
 */













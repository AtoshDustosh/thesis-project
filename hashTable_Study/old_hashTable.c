#include "old_hashTable.h"


static void _initHashTableTest();
static void _addHashCellTest();

static char* randomString(uint64_t strLength);
static uint64_t MyHash(char* str);
//static uint64_t BPHash(char* str);
//static uint64_t ELFhash(char *str);

void _HashTableTestSet() {
    _initHashTableTest();
    _addHashCellTest();
}

/*
 * Tests for working functions.
 */

/**
 * Test function initHashTable.
 */
static void _initHashTableTest() {
    printf("\n**************** _initHashTableTest ****************\n");
    printf("sizeof(HashCell): %d\n", sizeof(HashCell));
    uint64_t tableSize = 10;
    HashTable* hashTable = initHashTable(tableSize);

    displayHashTable(hashTable, tableSize);
}

/**
 * Test function addHashCell.
 */
static void _addHashCellTest() {
    printf("\n**************** _addHashCellTest ****************\n");
    uint64_t tableSize = 500000;
    HashTable* hashTable = initHashTable(tableSize);
    printf("hash table size: %"PRIu64"\n", tableSize);

    uint64_t i = 0;
    char* str;
    uint64_t strLength = 20;
    for(i = 0; i < tableSize; i++) {
        str = randomString(strLength);
//        printf("%"PRIu64": %s -> hash -> %"PRIu64"\n", i, str, MyHash(str));
        addHashCell(str, i, hashTable, tableSize);
        free(str);
    }

    printf("search hash cell of \"aaccgaca\" -> hash index: %"PRIu64"\n",
           searchHashIndexOfString("aaccgaca", hashTable, tableSize));
//    displayHashTable(hashTable, tableSize);
    checkHashTablePerformance(hashTable, tableSize);
}





/*
 * Working functions.
 */

HashCell* searchHashIndexOfString(char* str, HashTable* hashTable, uint64_t tableSize) {
    return hashTable->hashList[hashIndex(str, tableSize)];
}


void checkHashTablePerformance(HashTable* hashTable, uint64_t tableSize) {
    uint64_t i = 0;

    uint64_t totalBucketLength = 0;
    float bucketUsedCount = 0;
    for(i = 0; i < tableSize; i++) {
        HashCell* hashCell = hashTable->hashList[i];
        if(hashCell != NULL){
            totalBucketLength = totalBucketLength + hashCell->flag;
            bucketUsedCount = bucketUsedCount + (hashCell->flag && 1);
        }
    }

    printf("total bucket length: %"PRIu64"\n", totalBucketLength);
    printf("average bucket length: %f\n",
           bucketUsedCount == 0 ? 0 : totalBucketLength / bucketUsedCount);
    printf("bucket usage ratio: %f\n", bucketUsedCount / tableSize);
}

void displayHashTable(HashTable* hashTable, uint64_t tableSize) {
    uint64_t i = 0;

    uint64_t maxBucketLength = 0;
    for(i = 0; i < tableSize; i++) {
        uint64_t tempLength = 0;
        printf("hash cell %"PRIu64"\t", i);
        HashCell* hashCell = hashTable->hashList[i];
        while(hashCell != NULL){
            printf("(%"PRIu64", %d, 0x%p)\t", hashCell->data, hashCell->flag,
                   hashCell->nextCell);
            tempLength++;
            if(hashCell->nextCell != NULL){
                printf("->\t");
                hashCell = hashCell->nextCell;
            } else {
                break;
            }
        }
        printf("\n");

        if(tempLength > maxBucketLength) {
            maxBucketLength = tempLength;
        }
    }
    printf("Maximum bucket length: %"PRIu64"\n", maxBucketLength);
}

void addHashCell(char* str, uint64_t data, HashTable* hashTable, uint64_t tableSize) {
    uint64_t index = hashIndex(str, tableSize);
    HashCell* hashCell = hashTable->hashList[index];
    HashCell* newHashCell = (HashCell*)malloc(sizeof(HashCell));
    if(newHashCell == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    newHashCell->data = data;
    newHashCell->flag = 1;
    newHashCell->nextCell = NULL;
    if(hashCell == NULL){
        hashTable->hashList[index] = newHashCell;
    } else {
        while(hashCell != NULL){
            hashCell->flag++;
            if(hashCell->nextCell != NULL){
                hashCell = hashCell->nextCell;
            }else {
                break;
            }
        }
        hashCell->nextCell = newHashCell;
    }
}

HashTable* initHashTable(uint64_t tableSize) {
    HashTable* hashTable = (HashTable*)malloc(sizeof(HashTable));
    if(hashTable == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    hashTable->hashList = (HashCell**)malloc(sizeof(HashCell*) * tableSize);
    if(hashTable->hashList == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    uint64_t i = 0;
    if(hashTable->hashList == NULL) {
        printf("System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < tableSize; i++) {
        hashTable->hashList[i] = NULL;
    }
    hashTable->tableSize = tableSize;
    return hashTable;
}

uint64_t hashIndex(char *str, uint64_t tableSize) {
    return MyHash(str) % tableSize;
}





/*
 * Static functions. (file-localized functions)
 */




/**
 * My hash function.
 *
 * total bucket length: 5000000
 * average bucket length: 1.590123
 * bucket usage ratio: 0.628882
 *
 * @param str input string
 * @param hash value
 */
static uint64_t MyHash(char* str) {
    uint64_t hash = 0;
    uint64_t seed = 131;

    while(*str != '\0') {
//        uint64_t x = hash & 0xf000000000000000;
//        hash = (hash << 2) + *str;
//        if(x != 0) {
//            hash = hash ^ (x >> 60);
//            hash = hash & ~x;
//        }
        hash = hash * seed + (*str);
        str++;  // use pointers to go forward on a string
    }
    return (hash & 0x7fffffffffffffff);
}

///**
// * BPHash function. \note aborted
// *
// * total bucket length: 500000
// * average bucket length: 8.020018
// * bucket usage ratio: 0.124688
// *
// * @param str input string
// * @param hash value
// */
//static uint64_t BPHash(char* str) {
//    uint64_t hash = 0;
//    while(*str != '\0') {
//        hash = (hash << 7) ^ (*str);
//        str++;
//    }
//    return hash;
//}
//
///**
// * ELFHash function.
// *
// * total bucket length: 5000000
// * average bucket length: 1.801436
// * bucket usage ratio: 0.555113
// *
// * @param str input string
// * @param hash value
// */
//static uint64_t ELFhash(char *str) {
//    uint64_t hash = 0;
//
//    while(*str != '\0') {
//        uint64_t x = hash & 0xf000000000000000;
//        hash = (hash << 4) + *str;
//        if(x != 0) {
//            hash = hash ^ (x >> 60);
//            hash = hash & ~x;
//        }
//        str++;  // use pointers to go forward on a string
//    }
//    return (hash & 0x7fffffffffffffff);
//}

/**
 * Generate a random string at a specific length.
 *
 * @param strLength length of string
 * @return str generated string
 */
static char* randomString(uint64_t strLength) {
    char* str = (char*)malloc(sizeof(char) * (strLength + 1));
    if(str == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    uint64_t i = 0;
    for(i = 0; i < strLength; i++) {
        char randChar = (rand() % 26) % 4 + 'a';
        str[i] = randChar;
    }
    str[i] = '\0';
    return str;
}



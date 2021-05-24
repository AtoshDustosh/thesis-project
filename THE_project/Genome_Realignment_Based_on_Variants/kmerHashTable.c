#include "kmerHashTable.h"

/*********************************************************************
 *                         Structures and Accessors
 ********************************************************************/

struct KmerHashCell {
  char *kmerString;
  int32_t kmerPos;
  KmerHashCell *next;  // pointer to the next cell who has the same hash value
};

struct KmerHashTable {
  int32_t kmerLength;
  int32_t tableSize;
  int8_t *lineSize;         // array of cells' numbers on a hash line
  KmerHashCell **cellList;  // pointers to hash cells
};

struct Iterator_Kmerhash {
  KmerHashTable *table;
  int32_t cellIdx;  // index of cell
  int32_t lineIdx;  // index of the record on a hash line
};

inline char *kmerHashCell_string(KmerHashCell *cell) {
  return strdup(cell->kmerString);
}

inline int32_t kmerHashCell_pos(KmerHashCell *cell) { return cell->kmerPos; }

inline KmerHashCell *kmerHashCell_next(KmerHashCell *cell) {
  return cell->next;
}

/*********************************************************************
 *                  Auxiliary Structures and Methods
 ********************************************************************/

/**
 * @brief Generate a random kmer of specific length.
 */
static char *randomKmer(int32_t length) {
  char *kmer = (char *)calloc(length + 1, sizeof(char));
  if (kmer == NULL) {
    printf("Error: System memory not enough. \n");
    exit(EXIT_FAILURE);
  }
  int32_t i = 0;
  for (i = 0; i < length; i++) {
    char randChar = rand() % 4 + 'a';
    kmer[i] = randChar;
  }
  kmer[i] = '\0';
  return kmer;
}

static inline int32_t hashCode(char *kmer, int32_t pos, KmerHashTable *table) {
  // TODO try to improve its performance
  static const int32_t hashSeed = 131;
  int32_t value = 1;
  while (*kmer != '\0') {
    value = value * hashSeed + (*kmer);
    kmer++;
  }
  value = value + pos;
  return value & (INT32_MAX - 1);
}

static int32_t hashIndex(char *kmer, int32_t pos, KmerHashTable *table) {
  return hashCode(kmer, pos, table) % table->tableSize;
}

static inline KmerHashCell *init_KmerHashCell(char *kmer, int32_t pos) {
  KmerHashCell *cell = (KmerHashCell *)malloc(sizeof(KmerHashCell));
  if (cell == NULL) {
    fprintf(stderr, "Error: no enough meory for a new hash cell.\n");
    exit(EXIT_FAILURE);
  }
  cell->kmerPos = pos;
  cell->kmerString = strdup(kmer);
  cell->next = NULL;
  return cell;
}

static inline bool kmerHashCell_equal(KmerHashCell *cell, char *kmer,
                                      int32_t pos) {
  bool condition1 = strcmp(kmer, cell->kmerString) == 0 ? true : false;
  bool condition2 = pos == cell->kmerPos ? true : false;
  return condition1 && condition2;
}

/**
 * @brief  Destroy a kmer-hash-cell and return its pointer to the next cell.
 * @retval Poiner to the cell's next cell
 */
static inline KmerHashCell *destroy_KmerHashCell(KmerHashCell *cell) {
  KmerHashCell *next = cell->next;
  free(cell->kmerString);
  free(cell);
  return next;
}

/*********************************************************************
 *                            Basic Functions
 ********************************************************************/

Iterator_Kmerhash *init_iterator_Kmerhash(KmerHashTable *table) {
  Iterator_Kmerhash *it =
      (Iterator_Kmerhash *)malloc(sizeof(Iterator_Kmerhash));
  it->cellIdx = -1;
  it->lineIdx = -1;
  it->table = table;
  return it;
}

char *iterator_Kermhash_string(Iterator_Kmerhash *it) {
  if (it->cellIdx < it->table->tableSize && it->cellIdx >= 0) {
    KmerHashCell *cell = it->table->cellList[it->cellIdx];
    int lineSize = it->table->lineSize[it->cellIdx];
    if (it->lineIdx < lineSize) {
      int lineIdx = 0;
      while (lineIdx != it->lineIdx) {
        cell = cell->next;
        lineIdx++;
      }
      return strdup(cell->kmerString);
    } else {
      return NULL;
    }
  } else {
    return NULL;
  }
}

int32_t iterator_Kmerhash_pos(Iterator_Kmerhash *it) {
  if (it->cellIdx < it->table->tableSize) {
    KmerHashCell *cell = it->table->cellList[it->cellIdx];
    int lineSize = it->table->lineSize[it->cellIdx];
    if (it->lineIdx < lineSize) {
      int lineIdx = 0;
      while (lineIdx != it->lineIdx) {
        cell = cell->next;
        lineIdx++;
      }
      return cell->kmerPos;
    } else {
      return -1;
    }
  } else {
    return -1;
  }
}

bool iterator_Kmerhash_next(Iterator_Kmerhash *it) {
  if (it->cellIdx < it->table->tableSize) {
    if (it->cellIdx < 0 && it->lineIdx < 0) {
      it->cellIdx = 0;
      it->lineIdx = -1;
    }
    KmerHashCell *cell = it->table->cellList[it->cellIdx];
    int lineSize = it->table->lineSize[it->cellIdx];
    if (it->lineIdx + 1 < lineSize) {  // temporary line not finished
      it->lineIdx = it->lineIdx + 1;
      return true;
    } else {  // temporary line finished
      int cellIdx = it->cellIdx + 1;
      cell = it->table->cellList[cellIdx];
      while (cell == NULL && cellIdx < it->table->tableSize) {
        cellIdx++;
        cell = it->table->cellList[cellIdx];
      }
      if (cell == NULL) {  // all records have been iterated
        return false;
      } else {  // found another un-iterated record
        it->cellIdx = cellIdx;
        it->lineIdx = 0;
        return true;
      }
    }
  } else {
    return false;
  }
}

void destroy_iterator_Kmerhash(Iterator_Kmerhash *it) { free(it); }

KmerHashTable *init_kmerHashTable(int32_t tableSize, int32_t kmerLength) {
  KmerHashTable *hashTable = (KmerHashTable *)malloc(sizeof(KmerHashTable));
  if (hashTable == NULL) {
    fprintf(stderr,
            "Error: no enough memory for a kmer-hash-table of size: %" PRId32
            "\n",
            tableSize);
    exit(EXIT_FAILURE);
  }

  hashTable->kmerLength = kmerLength;
  hashTable->tableSize = tableSize;
  hashTable->lineSize = (int8_t *)calloc(tableSize, sizeof(int8_t));
  hashTable->cellList =
      (KmerHashCell **)calloc(tableSize, sizeof(KmerHashCell));
  if (hashTable->cellList == NULL || hashTable->lineSize == NULL) {
    fprintf(stderr,
            "Error: no enough memory for a kmer-hash-table of size: %" PRId32
            "\n",
            tableSize);
    exit(EXIT_FAILURE);
  }
  return hashTable;
}

void destroy_kmerHashTable(KmerHashTable *table) {
  for (int i = 0; i < table->tableSize; i++) {
    KmerHashCell *cell = table->cellList[i];
    while (cell != NULL) {
      cell = destroy_KmerHashCell(cell);
    }
  }
  free(table->lineSize);
  free(table);
}

bool kmerHashTable_add(char *kmer, int32_t pos, KmerHashTable *table) {
  if (strlen(kmer) != table->kmerLength) {
    fprintf(stderr,
            "Warning: incompatible length of kmer %ld. Should be %" PRId32
            ". \n",
            strlen(kmer), table->kmerLength);
    return false;
  }
  int32_t index = hashIndex(kmer, pos, table);
  KmerHashCell *cell = table->cellList[index];
  if (cell == NULL) {  // If this hash line is empty
    table->cellList[index] = init_KmerHashCell(kmer, pos);
  } else {  // If this hash line is not empty
    KmerHashCell *lastCell = NULL;
    while (cell != NULL && kmerHashCell_equal(cell, kmer, pos) == false) {
      lastCell = cell;
      cell = cell->next;
    }
    if (cell == NULL) {  // No identical cell found
      lastCell->next = init_KmerHashCell(kmer, pos);
    } else {         // Found identical cell
      return false;  // Ignore this kmer - duplicated
    }
  }
  table->lineSize[index]++;
  return true;
}

bool kmerHashTable_exist(char *kmer, int32_t pos, KmerHashTable *table) {
  int32_t index = hashIndex(kmer, pos, table);
  KmerHashCell *cell = table->cellList[index];
  if (cell == NULL) {  // If this hash line is empty
    return false;
  } else {  // If this hash line is not empty
    KmerHashCell *lastCell = NULL;
    while (cell != NULL && kmerHashCell_equal(cell, kmer, pos) == false) {
      lastCell = cell;
      cell = cell->next;
    }
    if (cell == NULL) {  // No identical cell found
      return false;
    } else {  // Found identical cell
      return true;
    }
  }
}

/*********************************************************************
 *                            Debug Functions
 ********************************************************************/

void print_kmerHashTable(KmerHashTable *table) {
  printf("table size: %" PRId32 ", kmer length: %" PRId32 "\n",
         table->tableSize, table->kmerLength);
  for (int32_t i = 0; i < table->tableSize; i++) {
    KmerHashCell *cell = table->cellList[i];
    printf("[%" PRId8 "]", table->lineSize[i]);
    while (cell != NULL) {
      printf("-> (%s,%" PRId32 ")", cell->kmerString, cell->kmerPos);
      cell = cell->next;
    }
    printf("\n");
  }
}

void statistics_kmerHashTable(KmerHashTable *table) {
  int32_t cnt_kmer = 0;
  int32_t cnt_occupied = 0;
  int max_length_occupied = 0;
  for (int32_t i = 0; i < table->tableSize; i++) {
    int8_t lineSize = table->lineSize[i];
    if (lineSize > 0) {
      cnt_occupied++;
      if (max_length_occupied < lineSize) max_length_occupied = lineSize;
    }
    cnt_kmer += lineSize;
  }

  printf("Count of stored kmers: %" PRId32 ".\n", cnt_kmer);
  printf("Occupied ratio: %f. \n", (float)cnt_occupied / table->tableSize);
  printf("Average length of occupied: %f. \n", (float)cnt_kmer / cnt_occupied);
  printf("Max length of occupied: %d. \n", max_length_occupied);
}

int _testSet_hashTable() {
  int32_t cnt_kmers = 200;
  int32_t kmerLength = 22;
  KmerHashTable *table = init_kmerHashTable(cnt_kmers, kmerLength);

  for (int32_t i = 0; i < cnt_kmers; i++) {
    char *kmer = randomKmer(kmerLength);
    int32_t pos = rand() & (INT32_MAX - 1);
    printf("kmer: %s, pos: %" PRId32 ", hashCode: %" PRId32 "\n", kmer, pos,
           hashCode(kmer, pos, table));
    kmerHashTable_add(kmer, pos, table);
    free(kmer);
  }

  Iterator_Kmerhash *it = init_iterator_Kmerhash(table);

  while (iterator_Kmerhash_next(it) != false) {
    printf("kmer: %s, pos: %d\n", iterator_Kermhash_string(it),
           iterator_Kmerhash_pos(it));
  }

  destroy_iterator_Kmerhash(it);

  print_kmerHashTable(table);

  statistics_kmerHashTable(table);

  destroy_kmerHashTable(table);
}
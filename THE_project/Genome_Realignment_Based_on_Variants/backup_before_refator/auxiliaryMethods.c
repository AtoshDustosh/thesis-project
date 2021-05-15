#include "auxiliaryMethods.h"

char *insertStr(const char *original, const char *inserted,
                const int insertPos) {
  int lengthOriginal = strlen(original);
  int lengthInserted = strlen(inserted);

  if (insertPos < 0 ||
      insertPos > lengthOriginal) {  // invalid argument "insertPos"
    return NULL;
  }

  // note that the "+1" is for the "\0" at the end of the string
  int lengthNewStr = lengthOriginal + lengthInserted;
  char *newStr = (char *)calloc(lengthNewStr + 1, sizeof(char));
  if (newStr == NULL) {
    fprintf(stderr, "Error: no enough memory for new string. \n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < insertPos; i++) {
    newStr[i] = original[i];
  }
  for (int i = 0; i < lengthInserted; i++) {
    newStr[i + insertPos] = inserted[i];
  }
  for (int i = 0; i < lengthOriginal - insertPos; i++) {
    newStr[insertPos + lengthInserted + i] = original[insertPos + i];
  }
  newStr[lengthInserted + lengthOriginal] = '\0';

  return newStr;
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/************************* Debug Methods ************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/

static int _test_insertStr() {
  const char *originals[] = {"012389", "56789", "012345", "123", "0123"};
  const char *inserted[] = {"4567", "01234", "6789", "0", "5"};
  const int pos[] = {4, 0, 6, -1, 5};
  const char *newStr[] = {"0123456789", "0123456789", "0123456789", NULL, NULL};

  int cnt = sizeof(newStr) / sizeof(char *);
  for (int i = 0; i < cnt; i++) {
    char *tmpStr = insertStr(originals[i], inserted[i], pos[i]);
    if (tmpStr == NULL) {
      assert(newStr[i] == NULL);
    } else {
      assert(strcmp(tmpStr, newStr[i]) == 0);
    }
    // printf("original: %s, inserted: %s\n", originals[i], inserted[i]);
    // printf("newStr: %s\n", tmpStr);
    free(tmpStr);
  }
  return 1;
}

void _testSet_auxiliaryMethods() { _test_insertStr(); }
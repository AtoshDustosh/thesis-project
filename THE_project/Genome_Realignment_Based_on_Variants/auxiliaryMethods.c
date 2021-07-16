#include "auxiliaryMethods.h"

char *revStr(const char *original, const int length) {
  if (original == NULL) {
    return NULL;
  }
  char *ret = (char *)calloc(length + 1, sizeof(char));

  for (int i = 0; i < length; i++) {
    ret[i] = original[length - i - 1];
  }
  ret[length] = '\0';

  return ret;
}

char *subStr(const char *original, const int begin, const int end) {
  if (original == NULL) {
    return NULL;
  }
  if (end < begin || begin < 0) {  // Invalid arguments
    return NULL;
  }
  char *ret = (char *)calloc((end - begin + 2), sizeof(char));
  if (ret == NULL) {
    fprintf(stderr, "Error: no enough memory for new string. \n");
    exit(EXIT_FAILURE);
  }

  int idx_original = 0;
  int idx_ret = 0;
  while (original[idx_original] != '\0') {
    if (idx_original >= begin && idx_original <= end) {
      ret[idx_ret++] = original[idx_original];
    }
    idx_original++;
  }
  ret[idx_ret] = '\0';

  return ret;
}

char *subStr_fast(const char *original, const int len_original, const int begin,
                  const int end) {
  if (original == NULL) {
    return NULL;
  }
  if (end < begin || begin < 0 || len_original <= 0 || end >= len_original) {
    return NULL;
  }

  char *ret = (char *)calloc((end - begin + 2), sizeof(char));
  if (ret == NULL) {
    fprintf(stderr, "Error: no enough memory for new string. \n");
    exit(EXIT_FAILURE);
  }

  int idx_original = begin;
  int idx_ret = 0;
  while (idx_original <= end) {
    ret[idx_ret++] = original[idx_original++];
  }
  ret[idx_ret] = '\0';

  return ret;
}

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

static int _test_revStr() {
  const char *originals[] = {"012389", "56789", "012345", "123", "0123"};
  int cnt = sizeof(originals) / sizeof(char *);
  int lengths[cnt];
  const char *reversed[] = {"983210", "98765", "543210", "321", "3210"};

  for (int i = 0; i < cnt; i++) {
    lengths[i] = strlen(originals[i]);
  }

  for (int i = 0; i < cnt; i++) {
    char *tmp_str = revStr(originals[i], lengths[i]);
    assert(tmp_str != NULL && strcmp(tmp_str, reversed[i]) == 0);
    free(tmp_str);
  }

  return 1;
}

static int _test_subStr() {
  const char *originals[] = {"012345", "012345", "012345",
                             "012345", "012345", "012345"};
  const int begins[] = {0, 0, 5, 2, 4, 3};
  const int ends[] = {5, 0, 5, 4, 3, 3};
  const char *subStrs[] = {"012345", "0", "5", "234", NULL, "3"};

  int cnt = sizeof(subStrs) / sizeof(char *);
  for (int i = 0; i < cnt; i++) {
    char *tmp_str = subStr(originals[i], begins[i], ends[i]);
    // printf("begin: %d, end: %d\n", begins[i], ends[i]);
    // printf("tmp str: %s, should str: %s\n", tmp_str, subStrs[i]);
    if (tmp_str == NULL) {
      assert(subStrs[i] == NULL);
    } else {
      assert(strcmp(subStrs[i], tmp_str) == 0);
    }
    free(tmp_str);
  }

  return 1;
}

static int _test_subStr_fast() {
  const char *originals[] = {"012345", "012345", "012345",
                             "012345", "012345", "012345"};
  const int len_originals[] = {6, 6, 6, 6, 6, 6};
  const int begins[] = {0, 0, 5, 2, 4, 3};
  const int ends[] = {5, 0, 5, 4, 3, 3};
  const char *subStrs[] = {"012345", "0", "5", "234", NULL, "3"};

  int cnt = sizeof(subStrs) / sizeof(char *);
  for (int i = 0; i < cnt; i++) {
    char *tmp_str = subStr_fast(originals[i], len_originals[i], begins[i], ends[i]);
    // printf("begin: %d, end: %d\n", begins[i], ends[i]);
    // printf("tmp str: %s, should str: %s\n", tmp_str, subStrs[i]);
    if (tmp_str == NULL) {
      assert(subStrs[i] == NULL);
    } else {
      assert(strcmp(subStrs[i], tmp_str) == 0);
    }
    free(tmp_str);
  }

  return 1;
}

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

void _testSet_auxiliaryMethods() {
  assert(_test_revStr());
  assert(_test_subStr());
  assert(_test_subStr_fast());
  assert(_test_insertStr());
}

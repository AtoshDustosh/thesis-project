#include "auxiliaryMethods.h"

char *insertStr(char *original, char *inserted, int insertPos) {
  int lengthOriginal = strlen(original);
  int lengthInserted = strlen(inserted);
  
  if(insertPos < 0 || insertPos > lengthOriginal){ // invalid argument "insertPos"
    return NULL;
  }
  
  // note that the "+1" is for the "\0" at the end of the string
  int lengthNewStr = lengthOriginal + lengthInserted;
  char *newStr = (char *)calloc(lengthNewStr + 1, sizeof(char));
  if (newStr == NULL) {
    fprintf(stderr, "Error: no enough memory for new string. \n");
    exit(EXIT_FAILURE);
  }

  int flagInserted = 0; // mark whether the insertion has been finished
  int idxOriginal = 0;
  int idxInserted = 0;
  for(int i = 0; i < lengthNewStr; i++){
    if(flagInserted == 0){
      // TODO
    }else{
      newStr[i] = original[idxOriginal++];
    }
  }
  

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
      assert(strcmp(tmpStr, newStr) == 0);
    }
    free(tmpStr);
  }
  return 1;
}

void _testSet_auxiliaryMethods() { _test_insertStr(); }
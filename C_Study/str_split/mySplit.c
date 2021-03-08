#include "mySplit.h"

static void _test_mySplit() {
  char *str = ">abs - def.gh-i -> j";
  char *delim1 = "-";
  char *delim2 = " - ";
  char *delim3 = ">";
  char *tokens1[] = {">abs ", " def.gh", "i ", "> j"};
  char *tokens2[] = {">abs", "def.gh-i -> j"};
  char *tokens3[] = {"abs - def.gh-i -", " j"};

  char **tokens = NULL;
  int tokenCnt = 0;

  tokenCnt = mySplit(str, delim1, &tokens);
  assert(tokenCnt == sizeof(tokens1) / sizeof(char *));
  for (int i = 0; i < tokenCnt; i++) {
    assert((strcmp(tokens[i], tokens1[i]) == 0) ||
           (printf("\"%s\", \"%s\"\n", tokens[i], tokens1[i]) >= 0));
  }

  tokenCnt = mySplit(str, delim1, &tokens);
  assert(tokenCnt == sizeof(tokens2) / sizeof(char *));
  for (int i = 0; i < tokenCnt; i++) {
    assert((strcmp(tokens[i], tokens2[i]) == 0) ||
           (printf("\"%s\", \"%s\"\n", tokens[i], tokens2[i]) >= 0));
  }

  tokenCnt = mySplit(str, delim2, &tokens);
  assert(tokenCnt == sizeof(tokens3) / sizeof(char *));
  for (int i = 0; i < tokenCnt; i++) {
    assert((strcmp(tokens[i], tokens3[i]) == 0) ||
           (printf("\"%s\", \"%s\"\n", tokens[i], tokens3[i]) >= 0));
  }
}

void _testSet_mySplit() { _test_mySplit(); }

int mySplit(char *str, char *delim, char ***tokens) {
  int tokenCnt = 0;
  // count tokens
  char *p_str = str;
  char *p_delim = delim;
  while (*p_str != '\0') {
    if (*p_str == *p_delim) {
      p_delim++;
    }
    if (*p_delim == '\0') {
      p_delim = delim;
      tokenCnt++;
    }
    p_str++;
  }

  // assign memory for tokens
  *tokens = (char **)malloc(sizeof(char *) * tokenCnt);
  if (*tokens == NULL) {
    fprintf(stderr, "Error: memory not enough. \n");
    exit(EXIT_FAILURE);
  }

  p_str = str;
  p_delim = delim;
  int tokenIdx = 0;
  int tokenSize = 0;
  while (*p_str != '\0') {
    if (*p_str == *p_delim) {
      p_delim++;
      tokenSize++;
    }
    if (*p_delim == '\0') {
      // this (tokenSize + 1) is in consideration of '\0'
      *tokens[tokenIdx] = (char *)malloc(sizeof(char) * (tokenSize + 1));
      if (*tokens[tokenIdx] == NULL) {
        fprintf(stderr, "Error: memory not enough. \n");
        exit(EXIT_FAILURE);
      }
      tokenIdx++;
      tokenSize = 0;
      p_delim = delim;
    }
    p_str++;
  }

  // copy tokens
  p_str = str;
  p_delim = delim;
  tokenIdx = 0;
  tokenSize = 0;
  // TODO
  while (*p_str != '\0') {
    if (*p_str == *p_delim) {
      p_delim++;
      tokenSize++;
    } else {
      *tokens[tokenIdx] = *p_str;
    }
    if (*p_delim == '\0') {
      tokenIdx++;
    }
    p_str++;
  }
  return 0;
}
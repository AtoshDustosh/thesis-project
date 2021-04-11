#include "autoPermutating.h"

Element *element_init(int recPos, int *recTypes, int *recLengths, int recCnt) {
  Element *e = (Element *)malloc(sizeof(Element));
  e->recCnt = recCnt;
  e->recPos = recPos;
  e->recTypes = (int *)calloc(recCnt, sizeof(int));
  e->recLengths = (int *)calloc(recCnt, sizeof(int));
  for (int i = 0; i < recCnt; i++) {
    e->recTypes[i] = recTypes[i];
    e->recLengths[i] = recLengths[i];
  }

  return e;
}

void element_print(Element *e) {
  for (int i = 0; i < e->recCnt; i++) {
    printf("[pos: %d, type: %d, length: %d]\n", e->recPos, e->recTypes[i],
           e->recLengths[i]);
  }
  printf("\n");
}

void element_destroy(Element *e) {
  free(e->recTypes);
  free(e->recLengths);
  free(e);
}

Combinations *combinations_init(int **combis, int combiSize, int combiCnt) {
  Combinations *cr = (Combinations *)malloc(sizeof(Combinations));
  cr->combis = combis;
  cr->combiSize = combiSize;
  cr->combiCnt = combiCnt;
  return cr;
}

void combinations_print(Combinations *cr) {
  if (cr == NULL) return;
  if (cr->combis == NULL) {
    return;
  }
  for (int i = 0; i < cr->combiCnt; i++) {
    for (int j = 0; j < cr->combiSize; j++) {
      printf("%d ", cr->combis[i][j]);
    }
    printf("\n");
  }
}

void combinations_destroy(Combinations *cr) {
  if (cr == NULL) return;
  if (cr->combiCnt == 0) {
    free(cr);
    return;
  }
  for (int i = 0; i < cr->combiCnt; i++) {
    free(cr->combis[i]);
  }
  free(cr->combis);
  free(cr);
}

EleCombinations *eleCombinations_init(int *eleCombi, int **alleleCombis,
                                      int combiSize, int combiCnt) {
  EleCombinations *ecr = (EleCombinations *)malloc(sizeof(EleCombinations));
  ecr->eleCombi = eleCombi;
  ecr->alleleCombis = alleleCombis;
  ecr->combiSize = combiSize;
  ecr->combiCnt = combiCnt;
  return ecr;
}

void eleCombinations_print(EleCombinations *ecr) {
  if (ecr == NULL) {
    fprintf(
        stderr,
        "Warning: null pointer occurred when printing an EleCOmbinations\n");
    return;
  }
  if (ecr->alleleCombis == NULL) {
    fprintf(
        stderr,
        "Warning: null pointer occurred when printing an EleCOmbinations\n");
    return;
  }
  for (int i = 0; i < ecr->combiCnt; i++) {
    printf("combis: ");
    for (int j = 0; j < ecr->combiSize; j++) {
      printf("(%d,%d) ", ecr->eleCombi[j], ecr->alleleCombis[i][j]);
    }
    printf("\n");
  }
}

void eleCombinations_destroy(EleCombinations *ecr) {
  if (ecr == NULL) return;
  if (ecr->alleleCombis == NULL) {
    free(ecr);
    return;
  }
  for (int i = 0; i < ecr->combiCnt; i++) {
    free(ecr->alleleCombis[i]);
  }
  free(ecr->alleleCombis);
  free(ecr);
}

static int recurseCombinations(int array[], int arraySize, int combi[],
                               int combiSize, int arrayIdx, int combiEleIdx,
                               int ***combis, int *combiIdx) {
  // combination with size 'combiSize' is selected
  if (combiEleIdx == combiSize) {
    if (*combis != NULL) {
      (*combis)[*combiIdx] = (int *)malloc(sizeof(int) * combiSize);
      for (int i = 0; i < combiSize; i++) {
        // printf("%d ", combi[i]);
        (*combis)[*combiIdx][i] = combi[i];
      }
      // printf("\n");
      *combiIdx = *combiIdx + 1;
    }
    return 1;
  }
  int combisCnt = 0;
  for (int i = arrayIdx; i < arraySize - (combiSize - combiEleIdx) + 1; i++) {
    combi[combiEleIdx] = array[i];
    combisCnt += recurseCombinations(array, arraySize, combi, combiSize, i + 1,
                                     combiEleIdx + 1, combis, combiIdx);
  }
  return combisCnt;
}

Combinations *combinations(int array[], int arraySize, int combiSize) {
  int **combis = NULL;
  int combiCnt = 0;

  int combiIdx = 0;
  int arrayIdx = 0;
  int combiEleIdx = 0;
  int *combi = (int *)calloc(combiSize, sizeof(int));
  combiCnt = recurseCombinations(array, arraySize, combi, combiSize, arrayIdx,
                                 combiEleIdx, &combis, &combiIdx);
  // printf("count of combinations: %d\n", combiCnt);
  combis = (int **)calloc(combiCnt, sizeof(int *));
  combiCnt = recurseCombinations(array, arraySize, combi, combiSize, arrayIdx,
                                 combiEleIdx, &combis, &combiIdx);

  free(combi);
  return combinations_init(combis, combiSize, combiCnt);
}

static void recurseEleCombinations(int *eleCombi, int combiSize,
                                   Element *eles[], int newAlleleIdx,
                                   int newAlleleCombi[], int **alleleCombis,
                                   int *eleRecCombiCnt) {
  if (newAlleleIdx == combiSize) {
    if (alleleCombis != NULL) {
      alleleCombis[*eleRecCombiCnt] = (int*)calloc(combiSize, sizeof(int));
      for(int i = 0; i < combiSize; i++){
        alleleCombis[*eleRecCombiCnt][i] = newAlleleCombi[i];
      }
    }
    // for (int i = 0; i < combiSize; i++) {
    //   printf("%d ", newAlleleCombi[i]);
    // }
    // printf("\n");
    *eleRecCombiCnt = *eleRecCombiCnt + 1;
    return;
  }
  // Iterate alleles of temporary elements which is under checking for new
  // allele
  for (int i = 0; i < eles[eleCombi[newAlleleIdx]]->recCnt; i++) {
    // If no element's allele has been selected, continue recursion
    if (newAlleleIdx == 0) {
      newAlleleCombi[newAlleleIdx] = i;
      recurseEleCombinations(eleCombi, combiSize, eles, newAlleleIdx + 1,
                             newAlleleCombi, alleleCombis, eleRecCombiCnt);
    } else {
      // If there exists selected element's allele, check if the last selected
      // allele covers the temporarily under checking allele. Holy shit ...
      if (eles[eleCombi[newAlleleIdx - 1]]->recPos +
              eles[eleCombi[newAlleleIdx - 1]]
                  ->recLengths[newAlleleCombi[newAlleleIdx - 1]] >
          eles[eleCombi[newAlleleIdx]]->recPos) {
        // if the last selected allele covers this allele
        continue;
      } else {
        // if the last selected allele doesn't cover this allele
        newAlleleCombi[newAlleleIdx] = i;
        recurseEleCombinations(eleCombi, combiSize, eles, newAlleleIdx + 1,
                               newAlleleCombi, alleleCombis, eleRecCombiCnt);
      }
    }
  }
  return;
}

EleCombinations *eleCombinations(int *eleCombi, int combiSize,
                                 Element *eles[]) {
  int **alleleCombis = NULL;

  int eleRecCombiCnt = 0;
  int newAlleleIdx = 0;
  int *newAlleleCombi = (int *)calloc(combiSize, sizeof(int));
  recurseEleCombinations(eleCombi, combiSize, eles, newAlleleIdx,
                         newAlleleCombi, alleleCombis, &eleRecCombiCnt);
  // TODO memset alleleCombis
  printf("ele allele combi cnt: %d\n", eleRecCombiCnt);
  alleleCombis = (int **)calloc(eleRecCombiCnt, sizeof(int *));
  eleRecCombiCnt = 0;
  recurseEleCombinations(eleCombi, combiSize, eles, newAlleleIdx,
                         newAlleleCombi, alleleCombis, &eleRecCombiCnt);

  free(newAlleleCombi);
  return eleCombinations_init(eleCombi, alleleCombis, combiSize,
                              eleRecCombiCnt);
}

// *********************************************************************
// *********************************************************************
// *********************************************************************
// ******************************** tests ******************************
// *********************************************************************
// *********************************************************************
// *********************************************************************

static void _test_Combinations() {
  int array[] = {0, 1, 2, 3, 4, 5};
  int arraySize = sizeof(array) / sizeof(int);
  for (int i = 1; i <= 6; i++) {
    Combinations *cr = combinations(array, arraySize, i);
    // combinations_print(cr);
    combinations_destroy(cr);
  }
}

static void _test_ElementPermutate() {
  /*
   *  positions of elements' records:
   *	1.1- 2.1- 3.1- 4.1- 5.1- 6.1-
   *       2.2-           5.2-
   *                      5.3-----
   *  1.2---------------------
   */
  const static int eleCnt = 6;
  int poss[] = {1, 5, 9, 12, 14, 19};
  int cnts[] = {2, 2, 1, 1, 3, 1};
  int lengths1[] = {1, 15};    // ele[0].rec[1] covers ele[1, 2, 3, 4]
  int lengths2[] = {3, 1};     // ele[1] doesn't cover ele[2]
  int lengths3[] = {3};        // ele[2] doesn't cover ele[3]
  int lengths4[] = {1};        // ele[3] doesn't cover ele[4]
  int lengths5[] = {1, 2, 6};  // ele[4].rec[2] covers ele[5]
  int lengths6[] = {2};
  int *lengths[] = {lengths1, lengths2, lengths3, lengths4, lengths5, lengths6};
  // These typeX are actually a marker for different records
  int types1[] = {11, 12};
  int types2[] = {21, 22};
  int types3[] = {31};
  int types4[] = {41};
  int types5[] = {51, 52, 53};
  int types6[] = {61};
  int *types[] = {types1, types2, types3, types4, types5, types6};

  Element *elements[eleCnt];
  for (int i = 0; i < eleCnt; i++) {
    elements[i] = element_init(poss[i], types[i], lengths[i], cnts[i]);
  }
  int elementsIdxes[] = {0, 1, 2, 3, 4, 5};

  int combiSize = 2;

  // Get all combinations of size "combiSize" from elements[]
  Combinations *cr = combinations(elementsIdxes, eleCnt, combiSize);
  // combinations_print(cr);
  // Get the advanced combinations of elements. Select the sub-elements.
  for (int i = 0; i < cr->combiCnt; i++) {
    printf("eleCombi: ");
    for (int j = 0; j < combiSize; j++) {
      printf("%d ", cr->combis[i][j]);
    }
    printf("\n");
    EleCombinations *ecr =
        eleCombinations(cr->combis[i], cr->combiSize, elements);
    eleCombinations_print(ecr);
    printf("\n");
    eleCombinations_destroy(ecr);
  }
  combinations_destroy(cr);

  for (int i = 0; i < eleCnt; i++) {
    element_destroy(elements[i]);
  }

  return;
}

static void _test_DoubleDimensionArray() {
  int arrayCnt = 6;
  int **arrays = (int **)malloc(sizeof(int *) * arrayCnt);
  int arraySizes[] = {2, 1, 1, 1, 3, 1};
  for (int i = 0; i < arrayCnt; i++) {
    arrays[i] = (int *)malloc(sizeof(int) * arraySizes[i]);
    for (int j = 0; j < arraySizes[i]; j++) {
      arrays[i][j] = (i + 1) * 10 + (j + 1);
    }
  }
  for (int i = 0; i < arrayCnt; i++) {
    for (int j = 0; j < arraySizes[i]; j++) {
      // printf("%d ", arrays[i][j]);
    }
    // printf("\n");
  }
  for (int i = 0; i < arrayCnt; i++) {
    free(arrays[i]);
  }
  free(arrays);
}

void _testSet_autoPermutating() {
  _test_Combinations();
  _test_DoubleDimensionArray();
  _test_ElementPermutate();
}
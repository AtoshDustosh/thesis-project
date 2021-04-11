#include "combinationsOfVars.h"

Combinations *combinations_init(int **combis, int combiSize, int combiCnt) {
  Combinations *cbs = (Combinations *)malloc(sizeof(Combinations));
  cbs->combis = combis;
  cbs->combiCnt = combiCnt;
  cbs->combiSize = combiSize;
  return cbs;
}

void combinations_print(Combinations *cbs) {
  if (cbs == NULL || cbs->combis == NULL) {
    fprintf(stderr,
            "Warning: null pointer occurred when printing a Combinations.\n");
    return;
  }
  for (int i = 0; i < cbs->combiCnt; i++) {
    for (int j = 0; j < cbs->combiSize; j++) {
      printf("%d ", cbs->combis[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void combinations_destroy(Combinations *cbs) {
  free(cbs->combis);
  free(cbs);
}

AlleleCombinations *init_AlleleCombinations(int *rvCombi, int **alleleCombis,
                                            int combiSize, int combiCnt) {
  AlleleCombinations *acbs =
      (AlleleCombinations *)malloc(sizeof(AlleleCombinations));
  acbs->rvCombi = rvCombi;
  acbs->alleleCombis = alleleCombis;
  acbs->combiSize = combiSize;
  acbs->combiCnt = combiCnt;
  return acbs;
}

void print_AlleleCombinations(AlleleCombinations *acbs) {
  if (acbs == NULL || acbs->rvCombi == NULL || acbs->alleleCombis == NULL) {
    fprintf(stderr,
            "Warning: null pointer occurred when printing an "
            "AlleleCombinations.\n");
    return;
  }
  printf("| rvCombi: ");
  for (int i = 0; i < acbs->combiSize; i++) {
    printf("%d ", acbs->rvCombi[i]);
  }
  printf("| alleleCombis: \n");
  for (int i = 0; i < acbs->combiCnt; i++) {
    printf("|\t");
    for (int j = 0; j < acbs->combiSize; j++) {
      printf("(%d,%d) ", acbs->rvCombi[j], acbs->alleleCombis[i][j]);
    }
    printf("\n");
  }
}

void destroy_AlleleCombinations(AlleleCombinations *acbs) {
  free(acbs->rvCombi);
  free(acbs->alleleCombis);
  free(acbs);
}

static int recurseCombinations(int array[], int arraySize, int newCombi[],
                               int combiSize, int arrayIdx, int combiEleIdx,
                               int **combis, int *combiIdx) {
  // combination with size 'combiSize' is selected
  if (combiEleIdx == combiSize) {
    if (combis != NULL) {
      combis[*combiIdx] = (int *)malloc(sizeof(int) * combiSize);
      for (int i = 0; i < combiSize; i++) {
        // printf("%d ", combi[i]);
        combis[*combiIdx][i] = newCombi[i];
      }
      // printf("\n");
      *combiIdx = *combiIdx + 1;
    }
    return 1;
  }
  int combisCnt = 0;
  for (int i = arrayIdx; i < arraySize - (combiSize - combiEleIdx) + 1; i++) {
    newCombi[combiEleIdx] = array[i];
    combisCnt += recurseCombinations(array, arraySize, newCombi, combiSize,
                                     i + 1, combiEleIdx + 1, combis, combiIdx);
  }
  return combisCnt;
}

Combinations *combinations(int array[], int arraySize, int combiSize) {
  int **combis = NULL;
  int combiCnt = 0;

  int combiIdx = 0;
  int arrayIdx = 0;
  int combiEleIdx = 0;
  int *newCombi = (int *)calloc(combiSize, sizeof(int));
  combiCnt = recurseCombinations(array, arraySize, newCombi, combiSize,
                                 arrayIdx, combiEleIdx, combis, &combiIdx);
  // printf("count of combinations: %d\n", combiCnt);
  combis = (int **)calloc(combiCnt, sizeof(int *));
  combiCnt = recurseCombinations(array, arraySize, newCombi, combiSize,
                                 arrayIdx, combiEleIdx, combis, &combiIdx);

  free(newCombi);
  return combinations_init(combis, combiSize, combiCnt);
}

static void recurseAlleleCombinations(ElementRecVcf *rvArray[], int rvCombi[],
                                      int combiSize, int newAlleleIdx,
                                      int newAlleleCombi[], int **alleleCombis,
                                      int *combiCnt) {
  if (newAlleleIdx == combiSize) {
    if (alleleCombis != NULL) {
      alleleCombis[*combiCnt] = (int *)calloc(combiSize, sizeof(int));
      for (int i = 0; i < combiSize; i++) {
        alleleCombis[*combiCnt][i] = newAlleleCombi[i];
      }
    }
    *combiCnt = *combiCnt + 1;
    return;
  }
  // Iterate alleles of rvs and select the next allele
  // "i = 1": ignore REF allele.
  // TODO this area is very likely to be buggy
  for (int i = 1; i < rvArray[rvCombi[newAlleleIdx]]->alleleCnt; i++) {
    // If no vcf record's allele has been selected, continue recursion.
    if (newAlleleIdx == 0) {
      newAlleleCombi[newAlleleIdx] =
          rvArray[rvCombi[newAlleleIdx]]->alleleIdx[i];
      recurseAlleleCombinations(rvArray, rvCombi, combiSize, newAlleleIdx,
                                newAlleleCombi, alleleCombis, combiCnt);
    } else {
      // If there exists selected vcf record's allele, check if the last
      // selected allele covers the temporarily under-checking allele.
      ElementRecVcf *lastSelectedErv = rvArray[rvCombi[newAlleleIdx - 1]];
      ElementRecVcf *tmpSelectedErv = rvArray[rvCombi[newAlleleIdx]];
      if (rvDataPos(lastSelectedErv->rv) +
          rvDataAlleleLength(lastSelectedErv->rv,
                             newAlleleCombi[newAlleleIdx - 1] >
                                 rvDataPos(tmpSelectedErv->rv))) {
        // if the last selected allele covers this allele
        continue;
      } else {
        newAlleleCombi[newAlleleIdx] =
            rvArray[rvCombi[newAlleleIdx]]->alleleIdx[i];
        recurseAlleleCombinations(rvArray, rvCombi, combiSize, newAlleleIdx + 1,
                                  newAlleleCombi, alleleCombis, combiCnt);
      }
    }
  }
  return;
}

AlleleCombinations *alleleCombinations(ElementRecVcf *ervArray[], int ervCnt,
                                       int rvCombi[], int combiSize) {
  // check rvCombi
  for (int i = 0; i < ervCnt; i++) {
    assert((rvCombi[i] < ervCnt) ||
           (fprintf(stderr, "Error: rvCnt incompatible with rvCombi. \n") < 0));
  }
  int **alleleCombis = NULL;
  int combiCnt = 0;

  int newAlleleIdx = 0;
  int *newAlleleCombi = (int *)calloc(combiSize, sizeof(int));
  // memeset alleleCombis
  recurseAlleleCombinations(ervArray, rvCombi, combiSize, newAlleleIdx,
                            newAlleleCombi, alleleCombis, &combiCnt);
  alleleCombis = (int **)calloc(combiCnt, sizeof(int *));
  combiCnt = 0;
  recurseAlleleCombinations(ervArray, rvCombi, combiSize, newAlleleIdx,
                            newAlleleCombi, alleleCombis, &combiCnt);

  free(newAlleleCombi);
  return init_AlleleCombinations(rvCombi, alleleCombis, combiSize, combiCnt);
}

// *********************************************************************
// *********************************************************************
// *********************************************************************
// ******************************** tests ******************************
// *********************************************************************
// *********************************************************************
// *********************************************************************

void _testSet_combinationsOfVars();
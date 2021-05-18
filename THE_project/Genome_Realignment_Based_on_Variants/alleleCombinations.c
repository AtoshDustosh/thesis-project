#include "alleleCombinations.h"

Combinations *init_combinations(int **combis, int length, int cnt) {
  Combinations *cbs = (Combinations *)malloc(sizeof(Combinations));
  cbs->combis = combis;
  cbs->length = length;
  cbs->cnt = cnt;
  return cbs;
}

void print_combinations(Combinations *cbs) {
  if (cbs == NULL || cbs->combis) {
    fprintf(stderr,
            "Warning: null pointer occurred when printing combinations.\n");
    return;
  }
  for (int i = 0; i < cbs->cnt; i++) {
    for (int j = 0; j < cbs->length; j++) {
      printf("%d ", cbs->combis[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void destroy_combinations(Combinations *cbs) {
  free(cbs);
}

Combinations_alleles *init_combination_alleles(int *combi_rv,
                                               int **combis_allele, int length,
                                               int cnt) {
  Combinations_alleles *acbs =
      (Combinations_alleles *)malloc(sizeof(Combinations_alleles));
  acbs->combi_rv = combi_rv;
  acbs->combis_allele = combis_allele;
  acbs->length = length;
  acbs->cnt = cnt;
  return acbs;
}

void print_combinations_alleles(Combinations_alleles *acbs) {
  if (acbs == NULL || acbs->combi_rv == NULL || acbs->combis_allele == NULL) {
    fprintf(stderr,
            "Warning: null pointer occurred when printing a "
            "Combinations_alleles.\n");
    return;
  }
  printf("| combi_rv: ");
  for (int i = 0; i < acbs->length; i++) {
    printf("%d ", acbs->combi_rv[i]);
  }
  printf("\n");
  printf("| combi_alleles (idx_erv, idx_allele_in_rv): \n");
  for (int i = 0; i < acbs->cnt; i++) {
    printf("|\t");
    for (int j = 0; j < acbs->length; j++) {
      printf("(%d,%d) ", acbs->combi_rv[j], acbs->combis_allele[i][j]);
    }
    printf("\n");
  }
}

void destroy_combinations_alleles(Combinations_alleles *acbs) {
  free(acbs);
}

static void recurse_combinations(int array[], int length_array, int combi_new[],
                                 int length_combi, int idx_array,
                                 int idx_combi_new, int **combis,
                                 int *cnt_combi) {
  // Combination with length "length_combi" selected
  if (idx_combi_new == length_combi) {
    if (combis != NULL) {
      combis[*cnt_combi] = (int *)calloc(length_combi, sizeof(int));
      for (int i = 0; i < length_combi; i++) {
        combis[*cnt_combi][i] = combi_new[i];
      }
    }
    *cnt_combi = *cnt_combi + 1;
    return;
  }
  for (int i = idx_array; i < length_array - (length_combi - idx_combi_new) + 1;
       i++) {
    combi_new[idx_combi_new] = array[i];
    recurse_combinations(array, length_array, combi_new, length_combi, i + 1,
                         idx_combi_new + 1, combis, cnt_combi);
  }
  return;
}

Combinations *calculate_combinations(int array[], int length_array,
                                     int length_combi) {
  int **combis = NULL;
  int cnt_combi = 0;

  int idx_array = 0;
  int idx_combi_new = 0;
  int combi_new[length_combi];
  recurse_combinations(array, length_array, combi_new, length_combi, idx_array,
                       idx_combi_new, combis, &cnt_combi);
  combis = (int **)calloc(cnt_combi, sizeof(int *));
  cnt_combi = 0;
  recurse_combinations(array, length_array, combi_new, length_combi, idx_array,
                       idx_combi_new, combis, &cnt_combi);
  return init_combinations(combis, length_combi, cnt_combi);
}

static void recurse_combinations_alleles(Element_RecVcf *rvArray[],
                                         int rvCombi[], int length_combi,
                                         int idx_combi_new, int combi_new[],
                                         int **combis, int *cnt_combi) {
  if (idx_combi_new == length_combi) {
    if (combis != NULL) {
      combis[*cnt_combi] = (int *)calloc(length_combi, sizeof(int));
      for (int i = 0; i < length_combi; i++) {
        combis[*cnt_combi][i] = combi_new[i];
      }
    }
    *cnt_combi = *cnt_combi + 1;
    return;
  }
  // This area is very likely to be buggy, but it actually works well
  for (int i = 0; i < rvArray[rvCombi[idx_combi_new]]->alleleCnt; i++) {
    if (idx_combi_new == 0) {
      // If there is no allele selected yet, continue recursion
      combi_new[idx_combi_new] = rvArray[rvCombi[idx_combi_new]]->alleleIdx[i];
      recurse_combinations_alleles(rvArray, rvCombi, length_combi,
                                   idx_combi_new + 1, combi_new, combis,
                                   cnt_combi);
    } else {
      // If some alleles have been selected, check if the last selected allele
      // coverse the temporarily under-checking allele
      Element_RecVcf *lastSelectedErv = rvArray[rvCombi[idx_combi_new - 1]];
      Element_RecVcf *tmpSelectedErv = rvArray[rvCombi[idx_combi_new]];
      if (rv_pos(lastSelectedErv->rv) +
              rv_alleleCoverLength(lastSelectedErv->rv,
                                   combi_new[idx_combi_new - 1]) >
          rv_pos(tmpSelectedErv->rv)) {
        // If the last selected allele covers this allele, ignore it
        continue;
      } else {
        combi_new[idx_combi_new] =
            rvArray[rvCombi[idx_combi_new]]->alleleIdx[i];
        recurse_combinations_alleles(rvArray, rvCombi, length_combi,
                                     idx_combi_new + 1, combi_new, combis,
                                     cnt_combi);
      }
    }
  }
}

Combinations_alleles *calculate_combinations_alleles(Element_RecVcf *ervArray[],
                                                     int length_array,
                                                     int *combi_rv,
                                                     int length_combi) {
  // Check rvCombi
  for (int i = 0; i < length_combi; i++) {
    assert((combi_rv[i] < length_array) ||
           (fprintf(
                stderr,
                "Error: combi_rv[%d](%d) incompatible with length_array(%d).\n",
                i, combi_rv[i], length_array) < 0));
  }
  int **combis = NULL;
  int cnt_combi = 0;

  int idx_combi_new = 0;
  int *combi_new = (int *)calloc(length_combi, sizeof(int));
  // memeset combis
  recurse_combinations_alleles(ervArray, combi_rv, length_combi, idx_combi_new,
                               combi_new, combis, &cnt_combi);
  combis = (int **)calloc(cnt_combi, sizeof(int *));
  cnt_combi = 0;
  recurse_combinations_alleles(ervArray, combi_rv, length_combi, idx_combi_new,
                               combi_new, combis, &cnt_combi);

  free(combi_new);
  return init_combination_alleles(combi_rv, combis, length_combi, cnt_combi);
}

void _testSet_alleleCombinations();
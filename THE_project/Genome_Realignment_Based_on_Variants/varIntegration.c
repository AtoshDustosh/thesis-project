#include "varIntegration.h"

void _testSet_varIntegration() {
  // TODO
}

void printVarIntegration(VarIntegration *vi) {
  printf("Variants to be integrated: \n");
  // TODO
  VarTodoChrom *vtc = vi->vtcs;
  static const uint32_t shiftLinePerCnt = 5;
  for (int i = 0; i < vi->chromCnt; i++) {
    assert(vtc != NULL ||
           (fprintf(stderr,
                    "Error: null pointer exception for VarTodoChrom when "
                    "printing VarIntegration.\n") >= 0));
    printf("chrom name: %s, variants count: %" PRIu32 "\n", vtc->name,
           vtc->varCnt);
    VarTodo *vt = vtc->vts->next;
    uint32_t printedCnt = 1;
    while (vt != NULL) {
      printf("%" PRIu32 "\t", vt->varIdx);
      if (printedCnt % shiftLinePerCnt == 0) {
        printf("\n");
      }
      printedCnt++;
      vt = vt->next;
    }
    if (printedCnt % shiftLinePerCnt != 1) {
      printf("\n");
    }
    vtc = vtc->next;
  }
}

static VarTodo *init_VarTodo() {
  VarTodo *vt = (VarTodo *)malloc(sizeof(VarTodo));
  vt->varIdx = -1;
  vt->next = NULL;
  return vt;
}

static void destroy_VarTodo(VarTodo *vt) { free(vt); }

static VarTodoChrom *init_VarTodoChrom() {
  VarTodoChrom *vtc = (VarTodoChrom *)malloc(sizeof(VarTodoChrom));
  vtc->name = "";
  vtc->varCnt = 0;
  vtc->vts = init_VarTodo();
  vtc->next = NULL;
  return vtc;
}

static void destroy_VarTodoChrom(VarTodoChrom *vtc) {
  VarTodo *vtTmp = vtc->vts;
  VarTodo *vtNxt = NULL;
  while (vtTmp != NULL) {
    vtNxt = vtTmp->next;
    destroy_VarTodo(vtTmp);
    vtTmp = vtNxt;
  }
  free(vtc);
}

VarIntegration *init_VarIntegration() {
  VarIntegration *vi = (VarIntegration *)malloc(sizeof(VarIntegration));
  vi->gv = NULL;
  vi->chromCnt = 0;
  vi->vtcs = NULL;
  return vi;
}

void destroy_VarIntegration(VarIntegration *vi) {
  VarTodoChrom *vtcTmp = vi->vtcs;
  VarTodoChrom *vtcNxt = NULL;
  while (vtcTmp != NULL) {
    vtcNxt = vtcTmp;
    destroy_VarTodoChrom(vtcTmp);
    vtcTmp = vtcNxt;
  }
  free(vi);
}

void addVtToVtc(uint32_t varIdx, VarTodoChrom *vtc) {
  /*
   * // TODO Consider making the variable "lastVt" static. That way you can use
   * it as a cache. And check whether varIdx is bigger than lastVt->varIdx
   * before tranversing down the linked-list. If all inputs are
   * varIdx-ascending, this will save all your time for traversing and speed up
   * the program. Oh! And don't forget to make a static variable to record the
   * last vtc. It would be bad if you use the same position on different
   * chromsomes.
   */
  VarTodo *tmpVt = vtc->vts->next;
  VarTodo *lastVt = vtc->vts;

  if (tmpVt == NULL) {
    VarTodo *vt = init_VarTodo();
    vt->varIdx = varIdx;
    vtc->vts->next = vt;
    vtc->varCnt++;
    return;
  }
  while (tmpVt != NULL) {
    if (tmpVt->varIdx < varIdx) {
      lastVt = tmpVt;
      tmpVt = tmpVt->next;
    } else if (tmpVt->varIdx == varIdx) {
      fprintf(stderr, "Warning: trying to add duplicated vt into a vtc.\n");
      return;
    } else {
      VarTodo *vt = init_VarTodo();
      vt->varIdx = varIdx;
      lastVt->next = vt;
      vt->next = tmpVt;
      vtc->varCnt++;
      return;
    }
  }
  // if varIdx is larger than any vt->varIdx in vtc
  VarTodo *vt = init_VarTodo();
  vt->varIdx = varIdx;
  lastVt->next = vt;
  vtc->varCnt++;
}

void addVtcToVarInt(VarTodoChrom *vtc, VarIntegration *vi) {
  VarTodoChrom *tmpVtc = vi->vtcs;
  if (tmpVtc == NULL) {
    vi->vtcs = vtc;
    vi->chromCnt++;
    return;
  }
  while (tmpVtc->next != NULL) {
    if (strcmp(tmpVtc->name, vtc->name) ==
        0) {  // if there already exists the same vtc
      fprintf(stderr, "Warning: trying to add duplicated VarTodoChrom. \n");
      return;
    }
    tmpVtc = tmpVtc->next;
  }
  // if same vtc not found, add it to the end of the linked-list of VarTodoChrom
  tmpVtc->next = vtc;
  vi->chromCnt++;
}

VarTodo getVtFromVtc(uint32_t vtIdx, VarTodoChrom *vtc) {
  if (vtIdx >= vtc->varCnt) {
    fprintf(stderr,
            "Error: array out of bound exception - vtIdx: %" PRIu32
            ", varCnt: %" PRIu32 "\n",
            vtIdx, vtc->varCnt);
    return NULL;
  }
  VarTodo *tmpVt = vtc->vts->next;
  for (int i = 0; i < vtIdx; i++) {
    tmpVt = tmpVt->next;
  }
  return tmpVt;
}

VarTodoChrom getVtcFromVarInt(char *chromName, VarIntegration *vi) {
  VarTodoChrom *tmpVtc = vi->vtcs;
  while (tmpVtc != NULL) {
    if (strcmp(chromName, tmpVtc->name) == 0) {
      return tmpVtc;
    } else {
      tmpVtc = tmpVtc->next;
    }
  }
  // if never found vtc with the same name as chromName
  return;
}
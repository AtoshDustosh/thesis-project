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

static VarTodoChrom *init_VarTodoChrom(){
  VarTodoChrom *vtc = (VarTodoChrom*)malloc(sizeof(VarTodoChrom));
  vtc->name = "";
  vtc->varCnt = 0;
  vtc->vts = init_VarTodo();
  vtc->next = NULL;
  return vtc;
}

static void destroy_VarTodoChrom(VarTodoChrom *vtc){
  VarTodo *vtTmp = vtc->vts;
  VarTodo *vtNxt = NULL;
  while(vtTmp != NULL){
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
  while(vtcTmp != NULL){
    vtcNxt = vtcTmp;
    destroy_VarTodoChrom(vtcTmp);
    vtcTmp = vtcNxt;
  }
  free(vi);
}

void addVtToVtc(VarTodo *vt, VarTodoChrom *vtc);

void addVtcToVarInt(VarTodoChrom *vtc, VarIntegration *vi);

VarTodo getVtFromVtc(uint32_t vtIdx, VarTodoChrom *vtc);

VarTodoChrom getVtcFromVarInt(char *chromName, VarIntegration *vi);

void markVar(char *chromName, uint32_t varIdx);


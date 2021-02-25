#include "C_revise.h"


static void _revise_Struct();

static void _revise_Struct(){
    printf("\n************ revise usage of struct *************\n");
    MyTuple *tuple = (MyTuple*)malloc(sizeof(MyTuple));
    tuple->name = "abcfdsafdsafdsa";
    tuple->value = 3;

    printf("my struct: {%d, %s}\n", tuple->value, tuple->name);
}

void CReviseTestSet(){
  _revise_Struct();
  printf("\n");
}

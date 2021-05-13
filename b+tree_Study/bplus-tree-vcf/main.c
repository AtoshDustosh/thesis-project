#include <stdio.h>
#include <stdlib.h>


#include "genomeVcf_bPlus.h"
#include "debug.h"



int main(int argc, char *argv[]){
  GenomeVcf_bplus *gv = genomeVcf_bplus_loadFile("test.vcf",4, 4);

  destroy_GenomeVcf_bplus(gv);

  return 0;
}
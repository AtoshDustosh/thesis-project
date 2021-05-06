#include "genomeVcf_bPlus.h"

GenomeVcf_bPlus *init_GenomeVcf_bPlus();

void destroy_GenomeVcf_bPlus(GenomeVcf_bPlus *gv_bPlus);

/************************************
 * Methods for manipulating GenomeVcf
 ************************************/
void genomeVcf_bPlus_insertRec(GenomeVcf_bPlus *gv_bPlus,
                               RecVcf_bPlus *rv_bPlus);

void genomeVcf_bPlus_removeRec(GenomeVcf_bPlus *gv_bPlus,
                               RecVcf_bPlus *rv_bPlus);

void genomeVcf_bPlus_loadFile(GenomeVcf_bPlus *gv_bPlus, char *filePath);

void genomeVcf_bPlus_writeFile(GenomeVcf_bPlus *gv_bPlus, char *filePath);

/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/************************* Debug Methods ************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/
/****************************************************************/

/**********************************
 * Debugging Methods for GenomeVcf
 **********************************/

void genomeVcf_bPlus_printRec(GenomeVcf_bPlus *gv_bPlus,
                              RecVcf_bPlus *rv_bPlus);

static int _test_BasicFunctions(){
  GenomeVcf_bPlus *gv_bPlus = init_GenomeVcf_bPlus();


  destroy_GenomeVcf_bPlus(gv_bPlus);
  return 1;
}

void _testSet_genomeVcf_bPlus(){
  assert(_test_BasicFunctions());
}
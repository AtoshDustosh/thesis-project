/*
 * @LastEditors: Atosh Dustosh
 */
#include "variantsStruct.h"

BcfData *initBcfData()
{
  BcfData *bD = (BcfData *)malloc(sizeof(BcfData));
  if (bD == NULL)
  {
    fprintf(stderr, "BcfData initialization failure - no enough memory. \n");
    exit(EXIT_SUCCESS);
  }
  return bD;
}

void destroyBcfData(BcfData *bD)
{
  bcf_hdr_destroy(bD->header);
  for (int i = 0; i < bD->n_variants; i++)
  {
    bcf_destroy1(bD->records[i]);
  }
  free(bD);
}
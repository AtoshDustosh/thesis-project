/*
 * @LastEditors: Atosh Dustosh
 */
#include "helperFunc.h"

int getFileType(const char *filePath)
{
  /*
   * fna, fasta | fastq | sam, bam | vcf, bcf
   * anf, atsaf | qtsaf | mas, mab | fcv, fcb
   */
  char suffix[MAX_SUFFIX_LENGTH] = {'\0'};
  int len = strlen(filePath);
  const char *nameCh = filePath + len - 1;

  for (int i = 0; i < len & *nameCh != '.'; i++)
  {
    suffix[i] = *nameCh;
    nameCh--;
  }
  if (strcmp(suffix, "anf") || strcmp(suffix, "atsaf"))
  {
    return FILE_REFERENCE;
  }
  else if (strcmp(suffix, "qtsaf"))
  {
    return FILE_READS;
  }
  else if (strcmp(suffix, "mas") || strcmp(suffix, "mab"))
  {
    return FILE_ALIGNMENTS;
  }
  else if (strcmp(suffix, "fcv") || strcmp(suffix, "fcb"))
  {
    return FILE_VARIANTS;
  }
  else
  {
    return FILE_UNKNOWN_TYPE;
  }
}
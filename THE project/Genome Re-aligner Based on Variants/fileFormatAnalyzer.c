/*
 * @LastEditors: Atosh Dustosh
 */
#include "fileFormatAnalyzer.h"

int getFileType(const char *filePath)
{
  /*
   * fna, fasta | fastq | sam, bam | vcf, bcf
   * anf, atsaf | qtsaf | mas, mab | fcv, fcb
   */
  char revSuffix[MAX_SUFFIX_LENGTH] = {'\0'};
  int len = strlen(filePath);
  const char *nameCh = filePath + len - 1;

  for (int i = 0; i < len & *nameCh != '.'; i++)
  {
    revSuffix[i] = *nameCh;
    nameCh--;
  }
  if (strcmp(revSuffix, "anf") == 0 || strcmp(revSuffix, "atsaf") == 0 || strcmp(revSuffix, "af") == 0)
  {
    return FILE_REFERENCE;
  }
  else if (strcmp(revSuffix, "qtsaf") == 0)
  {
    return FILE_READS;
  }
  else if (strcmp(revSuffix, "mas") == 0 || strcmp(revSuffix, "mab") == 0)
  {
    return FILE_ALIGNMENTS;
  }
  else if (strcmp(revSuffix, "fcv") == 0 || strcmp(revSuffix, "fcb") == 0)
  {
    return FILE_VARIANTS;
  }
  else
  {
    return FILE_UNKNOWN_TYPE;
  }
}
/*
 * @Author: Atosh Dustosh
 * @Date: 2021-01-28 21:02:40
 * @LastEditTime: 2021-01-28 21:54:41
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /Genome Re-aligner Based on Variants/variantsStruct.h
 */
#ifndef VARIANTSSTRUCT_H_INCLUDED
#define VARIATNSSTRUCT_H_INCLUDED

#include <htslib/vcf.h>

/*
 * These 2 aliases are used only outside the header file and cannot be 
 * used when implementing header functions. They only exist for better
 * reading. 
 */
typedef bcf_hdr_t *BcfHeader;
typedef bcf1_t *BcfRecord;

typedef struct _define_BcfData
{
  int n_variants;
  BcfHeader header;
  BcfRecord *records;
} BcfData;

/**
 * @brief  Allocate and initialize a BcfData object. 
 * @note   The BcfData struct returned by a successful call should be
 * freed via destroyBcfData() when it is no longer needed. 
 * @retval pointer to an empty BcfData; NULL otherwise
 */
BcfData *initBcfData();

/**
 * @brief  Deallocate a BcfData object.
 */
void destroyBcfData(BcfData *bD);

#endif // VARIATNSSTRUCT_H_INCLUDED
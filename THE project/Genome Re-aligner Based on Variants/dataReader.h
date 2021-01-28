/*
 * @Author: your name
 * @Date: 2021-01-28 21:26:34
 * @LastEditTime: 2021-01-28 22:08:10
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /Genome Re-aligner Based on Variants/dataReader.h
 */
#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "convenience.h"
#include "variantsStruct.h"

/**
 * @brief  Read the bcf/vcf file and load data into a BcfData object. 
 * @param  *opts: Options got from console
 * @param  bD: bcf data collecion object
 * @retval None
 */
// TODO I'm not sure if I should use "Options" to pass the inputFilePath.
// But in case there are other requirements when reading data, I'll just
// leave it be.
void readBcfFile(const Options *opts, BcfData *bD);

#endif
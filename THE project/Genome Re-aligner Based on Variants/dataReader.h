/*
 * @Author: your name
 * @Date: 2021-01-28 21:26:34
 * @LastEditTime: 2021-01-29 20:05:31
 * @LastEditors: Atosh Dustosh
 * @Description: In User Settings Edit
 * @FilePath: /Genome Re-aligner Based on Variants/dataReader.h
 */
#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "convenience.h"
#include "fileFormatAnalyzer.h"
#include "variantsStruct.h"


/**
 * @brief  Glimpse a file. 
 * @note  This method is designed for checking the
 * contents of a very large file that cannot be opened directly. 
 */
void glimpseFile(const Options *opts);

/**
 * @brief  Read the bcf/vcf file and load data into a BcfData object. 
 */
// TODO I'm not sure if I should use "Options" to pass the inputFilePath.
// But in case there are other requirements when reading data, I'll just
// leave it be.
void readBcfFile(const Options *opts, BcfData *bD);

#endif
#ifndef GRBVOPERATION_H_INCLUDED
#define GRBVOPERATION_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <inttypes.h>

#include "alignment.h"
#include "auxiliaryMethods.h"
#include "debug.h"
#include "genomeFa.h"
#include "genomeSam.h"
#include "grbvOptions.h"

/**
 * @brief select reads with MAPQ lower than MAPQ_threshold from previously set
 * sam file and then output them into the previously set output file
 */
void selectBadReads(Options *opts);

/**
 * @brief  Collect statistics from the input vcf file and print into the
 * specified output file, or if not specified, to the console.
 */
void statistics_vcf(Options *opts);

void _testSet_grbvOperations();

#endif
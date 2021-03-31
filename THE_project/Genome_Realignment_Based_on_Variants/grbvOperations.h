#ifndef GRBVOPERATION_H_INCLUDED
#define GRBVOPERATION_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <inttypes.h>

#include "alignment.h"
#include "auxiliaryMethods.h"
#include "genomeFa.h"
#include "genomeSam.h"
#include "genomeVcf.h"
#include "grbvOptions.h"

/**
 * @brief  Integrate variants in *.vcf files into alignment records of reads in
 * *.sam files. Do realignment during the integration and modify the cigars and
 * mapqs of them. And output the integrated results into the designated outpt
 * file.
 */
void integrateVcfToSam(Options *opts);

/**
 * @brief  A refactored version of the method "integrateVcfToSam". The original
 * implementation is TOO redundant and stupid, and the result is not exactly
 * what we want.
 */
void integrateVcfToSam_refactored(Options *opts);

/**
 * @brief select reads with MAPQ lower than MAPQ_threshold from previously set
 * sam file and then output them into the previously set output file
 */
void selectBadReads(Options *opts);

#endif
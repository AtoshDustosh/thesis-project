#ifndef GRBVOPERATION_H_INCLUDED
#define GRBVOPERATION_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <inttypes.h>

#include "alignment.h"
#include "genomeFa.h"
#include "genomeSam.h"
#include "genomeVcf.h"
#include "grbvOptions.h"

/**
 * @brief select reads with MAPQ lower than MAPQ_threshold from previously set
 * sam file and then output them into the previously set output file
 */
void selectBadReads(Options *opts);

/**
 * @brief  Integrate variants in *.vcf files into alignment records of reads in
 * *.sam files. Do realignment during the integration and modify the cigars and
 * mapqs of them.
 */
void integrateVcfToSam(Options *opts);

#endif
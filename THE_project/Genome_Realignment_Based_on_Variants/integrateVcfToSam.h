#ifndef INTEGRATEVCFTOSAM_H_INCLUDED
#define INTEGRATEVCFTOSAM_H_INCLUDED

#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <inttypes.h>
#include <pthread.h>

#include "alignment.h"
#include "auxiliaryMethods.h"
#include "alleleCombinations.h"
#include "debug.h"
#include "genomeFa.h"
#include "genomeSam.h"
#include "genomeVcf_bPlus.h"
#include "grbvOptions.h"


static char* default_outputFile = "defaultOutput.txt";

/**
 * @brief  Final version of integrating vcf records into sam records.
 */
void integration(Options* opts);

#endif
#ifndef DATAREADER_H_INCLUDED
#define DATAREADER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "myMacros.h"
#include "infoPrinter.h"

/**
 * Scan the input file and print all sam records using console.
 *
 * @param inputFile
 * @param opts
 * @return status
 */
int sam_scan_and_print(htsFile *inputFile, Options *opts);

/**
 * Scan the input file and print all vcf records using console.
 *
 * @param inputFile
 * @param opts
 * @return status
 */
int vcf_scan_and_print(htsFile *inputFile, Options *opts);


#endif // DATAREADER_H_INCLUDED

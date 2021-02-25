#ifndef INFOPRINTER_H_INCLUDED
#define INFOPRINTER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>


/**
 * Print information about the header of a vcf/bcf file.
 */
void printVcfHeader(bcf_hdr_t *header);

/**
 * Print a vcf/bcf record on the console.
 */
void printVcfRecord(bcf1_t *record);

/**
 * Print information about the header of a sam/bam file.
 */
void printSamHeader(bam_hdr_t *header);

/**
 * Print a sam/bam record in detail on the console.
 */
void printSamRecord(bam1_t *record);


#endif // INFOPRINTER_H_INCLUDED

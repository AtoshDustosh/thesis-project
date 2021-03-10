#ifndef INFOPRINTER_H_INCLUDED
#define INFOPRINTER_H_INCLUDED

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Print information about the header of a vcf/bcf file.
 */
void printVcfHeader(bcf_hdr_t *header);

/**
 * @brief Print detailed information with annotations. Do not use this for
 * printing huge amount of vcf records, but use "printVcfRecord_brief" instead.
 */
void printVcfRecord(bcf1_t *record);

/**
 * @brief  Print vcf record using standard format of *.vcf files.
 */
void printVcfRecord_brief(bcf_hdr_t *hdr, bcf1_t *record);

/**
 * Print information about the header of a sam/bam file.
 */
void printSamHeader(bam_hdr_t *header);

/**
 * @brief Print detailed information with annotations. Do not use this for
 * printing huge amount of sam records, but use "printSamRecord_brief" instead.
 */
void printSamRecord(bam1_t *record);

/**
 * @brief  Print sam record using standard format of *.sam files.
 */
void printSamRecord_brief(bam_hdr_t *hdr,bam1_t *record);

#endif  // INFOPRINTER_H_INCLUDED

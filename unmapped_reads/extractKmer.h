#ifndef EXTRACTKMER_H_INCLUDED
#define EXTRACTKMER_H_INCLUDED

#pragma once

#include <htslib/vcf.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief  Extract kmers.
 * @param  *interval: sequence from reference genome
 * @param  start_interval: 1-based start position of "*interval"
 * @param  end_interval: 1-based end position of interval
 * @param  *variants[]: array of variants that need to be integrated
 * @param  cnt_variants: count of variants / length of "variants[]"
 * @param  length_kmer: length of a kmer
 * @retval  *ret_cnt_kmer: count of kmers found
 * @retval  *ret_kmers[]: kmers
 * @retval  *ret_pos_kmers[]: start positions of kmers
 */
void extractKmer(const char *interval, int64_t start_interval,
                 int64_t end_interval, bcf1_t *variants[], int cnt_variants,
                 int length_kmer, int *ret_cnt_kmer, char *ret_kmers[],
                 int64_t *ret_pos_kmers[]);

#endif
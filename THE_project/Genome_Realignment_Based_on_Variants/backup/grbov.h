/*
 * @LastEditors: Atosh Dustosh
 */
#ifndef GRBOV_H_INCLUDED
#define GRBOV_H_INCLUDED

#pragma once

/**
 * @brief  Select bad reads from input file *.sam/bam and then write them into output file *.sam/bam with necessary information. 
 */
void selectBadReads(const Options *opts);

#endif

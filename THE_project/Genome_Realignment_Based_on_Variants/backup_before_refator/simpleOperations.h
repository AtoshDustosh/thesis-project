#ifndef SIMPLEOPERATIONS_H_INCLUDED
#define SIMPLEOPERATIONS_H_INCLUDED

#pragma once

#include "grbvOptions.h"

/**
 * @brief Count the records in all designated files.
 */
void countRec(Options *opts);

/**
 * @brief Print the first [number] lines of all files to the console.
 */

void firstLines(Options *opts);

/**
 * @brief  Extract bases of the selected chromosome and write into designated
 * output file together with the chromosome's info field.
 */
void extractChrom(Options *opts);

#endif
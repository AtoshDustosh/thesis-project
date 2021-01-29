/* 
 * This header file must contain only simple defination and structures. 
 * No subsequent process on defination and structures allowed. 
 * @LastEditors: Atosh Dustosh
 */
#ifndef CONVENIENCE_H_INCLUDED
#define CONVENIENCE_H_INCLUDED

#pragma once

#include <string.h>

// File type
#define FILE_UNKNOWN_TYPE 0
#define FILE_VARIANTS 1
#define FILE_REFERENCE 2
#define FILE_ALIGNMENTS 3
#define FILE_READS 4

// Operation types of Usage.
#define NO_OPERATION 0
#define OPT_GLIMPSE 1

// Operation types of glimpse
#define GLIMPSE_NOT_DESIGNATED 0
#define GLIMPSE_FIRST_LINES 1
#define GLIMPSE_COUNT_RECORDS 2

typedef struct _define_GlimpseArgs
{
  int type;
  int n_firstLines;
} GlimpseArgs;

typedef struct _define_Options
{
  char *inputFilePath;
  char *outputFilePath;
  int operationType;
  GlimpseArgs glimpse;
} Options;


#endif // CONVENIENCE_H_INCLUDED

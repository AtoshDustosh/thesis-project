/* 
 * This header file must contain only simple defination and structures. 
 * No subsequent process on defination and structures allowed. 
 * @LastEditors: Atosh Dustosh
 */
#ifndef GRBVOPTIONS_H_INCLUDED
#define GRBVOPTIONS_H_INCLUDED

#pragma once

#include <stdio.h>
#include <stdlib.h>

/*
 * Operation types of Usage.
 */
#define NO_OPERATION 0

/*
 * Set input and output fies. 
 */
#define OPT_SET_OUTPUTFILE 1
#define OPT_SET_FAFILE 2
#define OPT_SET_FASTQFILE 3
#define OPT_SET_SAMFILE 4
#define OPT_SET_VCFFILE 5

/*
 * Some simple operations against files. 
 */
#define OPT_COUNTREC 6
#define OPT_FIRSTLINES 7

/*
 * GRBV operations.
 */
#define OPT_SELECTBADREADS 8
#define OPT_COMPARESAM 9

typedef struct _define_Options
{
  char *faFile;
  char *fastqFile;
  char *samFile;
  char *vcfFile;
  char *outputFile;

  int countRec;
  int firstLines;

  int selectBadReads;
  int compareSam;
} Options;

typedef struct _define_FileList
{
  char **paths;
  int count;
} FileList;

/**
 * @brief Get the files designated by command inputs. 
 * 
 * @param opts command inputs
 * @return FileList* a list of designated files. You can access the value "count" to get the size of it. The list must be freed mannually later with destroyFileList(). 
 */
FileList* designatedFiles(Options *opts);

/**
 * @brief Destroy the FileList object
 */
void destroyFileList(FileList* fl);

#endif // GRBVOPTIONS_H_INCLUDED

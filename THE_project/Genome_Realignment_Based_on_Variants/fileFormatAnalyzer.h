/*
 * @LastEditors: Atosh Dustosh
 */

#ifndef FILEFORMATANALYZER_H_INCLUDED
#define FILEFORMATANALYZER_H_INCLUDED

#pragma once

#include <string.h>


// File type
#define FILE_UNKNOWN_TYPE 0
#define FILE_VARIANTS 1
#define FILE_REFERENCE 2
#define FILE_ALIGNMENTS 3
#define FILE_READS 4

// maximum file format suffix name length
#define MAX_SUFFIX_LENGTH 6

#define MAX_LINE_BUFFER 2048


int getFileType(const char *filePath);

#endif // FILEFORMATANALYZER_H_INCLUDED
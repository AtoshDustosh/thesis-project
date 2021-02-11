/*
 * @Author: Atosh Dustosh
 * @Date: 2021-01-28 09:23:44
 * @LastEditTime: 2021-01-29 20:00:51
 * @LastEditors: Atosh Dustosh
 * @Description: In User Settings Edit
 * @FilePath: /Genome Re-aligner Based on Variants/main.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <htslib/vcf.h>

#include "convenience.h"
#include "dataReader.h"
#include "variantsStruct.h"

const char *optStr = "i:o:gF:C";
static struct option optInitArray[] = {
    {"output", required_argument, NULL, 'o'},
    {"input", required_argument, NULL, 'i'},
    {"glimpse", no_argument, NULL, 'g'},
    {"firstlines", required_argument, NULL, 'F'},
    {"count", no_argument, NULL, 'C'},
    {0, 0, 0, 0},
};

void Usage()
{
    printf("Usage: program <command> <arguments>\n");
    printf("\n");
    printf("Options:\n");
    printf("\t<-i|--input> [input file] set input file\n");
    printf("\t<-o|--output> [output file] set output file\n");
    printf("\t<-g|--glimpse> glimpse the input file\n");
    printf("\t\t<-F|--firstlines> [integer number]"
           "\t\tscan the first [integer number] lines of input file\n");
    printf("\t\t<-C|--count> count record number or lines of input file\n");
}

int main(int argc, char *argv[])
{
    Options options;

    options.inputFilePath = NULL;
    options.outputFilePath = NULL;
    options.operationType = NO_OPERATION;
    options.glimpse.type = GLIMPSE_NOT_DESIGNATED;
    options.glimpse.n_firstLines = 0;

    int optRet = 0;
    while ((optRet = getopt_long(argc, argv, optStr, optInitArray, NULL)) != -1)
    {
        switch (optRet)
        {
        // basic operations
        case 'i':
        {
            printf("input file: %s\n", optarg);
            options.inputFilePath = optarg;
            break;
        }
        case 'o':
        {
            printf("output file: %s\n", optarg);
            options.outputFilePath = optarg;
            break;
        }
        case 'g':
        {
            if (options.operationType != NO_OPERATION)
            {
                fprintf(stderr, "Error: multiple options conflicted.\n");
                return EXIT_SUCCESS;
            }
            else
            {
                options.operationType = OPT_GLIMPSE;
            }
            break;
        }
        // subsequent operations
        case 'F':
        {
            if (options.glimpse.type != GLIMPSE_NOT_DESIGNATED)
            {
                fprintf(stderr, "Error: multiple glimpse options conflicted. \n");
                return EXIT_SUCCESS;
            }
            else
            {
                options.glimpse.type = GLIMPSE_FIRST_LINES;
                options.glimpse.n_firstLines = atoi(optarg);
            }
            break;
        }
        case 'C':
        {
            if (options.glimpse.type != GLIMPSE_NOT_DESIGNATED)
            {
                fprintf(stderr, "Error: multiple glimpse options conflicted. \n");
                return EXIT_SUCCESS;
            }
            else
            {
                options.glimpse.type = GLIMPSE_COUNT_RECORDS;
            }
            break;
        }
            // TODO other options to be realized
        }
    }

    /*
     * The switch block above sets the basic and subsequent operations,
     * but does no process on the data files. The switch block below is
     * the actual processing block. 
     */
    switch (options.operationType)
    {
    case NO_OPERATION:
        printf("No operation type selected. \n");
        Usage();
        return EXIT_SUCCESS;
    case OPT_GLIMPSE:
        printf("Glipmse file %s\n", options.inputFilePath);
        glimpseFile(&options);
        break;
    default:
        printf("Unknown operation.\n");
        Usage();
        return EXIT_SUCCESS;
    }
    return EXIT_SUCCESS;
}
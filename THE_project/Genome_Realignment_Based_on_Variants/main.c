/*
 * @Author: Atosh Dustosh
 * @Date: 2021-01-28 09:23:44
 * @LastEditTime: 2021-02-21 21:07:33
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
#include "grbov.h"

const char *optStr = "i:o:gF:C";
int lopt = 0;
static struct option optInitArray[] = {
    {"output", required_argument, NULL, 'o'},
    {"input", required_argument, NULL, 'i'},
    {"glimpse", no_argument, NULL, 'g'},
    {"firstlines", required_argument, NULL, 'F'},
    {"count", no_argument, NULL, 'C'},
    {"selectbadreads", no_argument, &lopt, 1},
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
    printf("\t\t<-F|--firstlines> [integer number] scan the first [integer number] lines of input file\n");
    printf("\t\t<-C|--count> count record number or lines of input file\n");
    printf("\t<--selectbadreads> select bad mapping reads from input file (*.sam/bam) and put them into the output file (*.sam/bam) with necessary infos.\n");
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
        // short options
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
    // long options
    case 0:
    {
        switch (lopt)
        {
        case 1:
        {
            if (options.operationType != NO_OPERATION)
            {
                fprintf(stderr, "Error: multiple options conflicted.\n");
                return EXIT_SUCCESS;
            }
            else
            {
                options.operationType = OPT_SELECT_BAD_READS;
            }
            break;
        }
        default:
            exit(EXIT_SUCCESS);
            break;
        }
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
    case OPT_SELECT_BAD_READS:
        printf("Select bad mapping reads from %s, and then output them into %s with necessary information", options.inputFilePath, options.outputFilePath);
        break;
    default:
        printf("Unknown operation.\n");
        Usage();
        return EXIT_SUCCESS;
    }
    return EXIT_SUCCESS;
}
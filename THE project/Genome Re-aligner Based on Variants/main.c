/*
 * @Author: Atosh Dustosh
 * @Date: 2021-01-28 09:23:44
 * @LastEditTime: 2021-01-28 22:17:23
 * @LastEditors: Atosh Dustosh
 * @Description: In User Settings Edit
 * @FilePath: /Genome Re-aligner Based on Variants/main.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <htslib/vcf.h>

#include "convenience.h"
 
struct option longOpts[] = {
    {"Output", required_argument, NULL, 'o'},
    {"Input", required_argument, NULL, 'i'},
    // {"love", required_argument, NULL, 'l'},
    {0, 0, 0, 0},
};


void Usage()
{
    printf("Usage: program <command> [options]\n");
    printf("\n");
    printf("Options:\n");
    printf("\t<-i|--Input> [input file] <-o|--Output> [output file]\n");
}


int main(int argc, char *argv[])
{
    Options options;
    int optRet = 0;
    char *optStr = "io:";
    while ((optRet = getopt_long(argc, argv, optStr, longOpts, NULL)) != -1)
    {
        switch (optRet)
        {
        case 'i':
            printf("input file: %s\n", optarg);
            options.inputFilePath = optarg;
            break;
        case 'o':
            printf("output file: %s\n", optarg);
            options.outputFilePath = optarg;
            break;
        // TODO other options to be realized
        default:
            Usage();
            return EXIT_SUCCESS;
        }
    }
    return EXIT_SUCCESS;
}
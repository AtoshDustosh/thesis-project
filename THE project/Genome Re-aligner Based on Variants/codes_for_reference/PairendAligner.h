#ifndef PAIRENDALIGNER_H_INCLUDED
#define PAIRENDALIGNER_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>

#include "SNAP.h"
#include "Read.h"

typedef struct _define_PairendAligner {
    SNAP* snap;
    char* fastqFilePath1;
    FILE* fpointer1;
    Read* read1;
    Read* read1_reverseCompliment;

    char* fastqFilePath2;
    FILE* fpointer2;
    Read* read2;
    Read* read2_reverseCompliment;
} PairendAligner;

/**
 * Load an initialized SNAP structure into pair-end aligner.
 *
 * @param snap a SNAP structure
 * @param pairendAligner the pair-end aligner
 */
void loadSNAPStructure(SNAP* snap, PairendAligner* pairendAligner);

/**
 * Initialize the pair-end aligner with 2 fastq file paths.
 *
 * @param fastqFilePath1 file path of fastq file 1
 * @param fastqFilePath2 file path of fastq file 2
 * @param pairendAligner the pair-end aligner
 */
void initPairendAligner(char* fastqFilePath1, char* fastqFilePath2, PairendAligner* pairendAligner);






#endif // PAIRENDALIGNER_H_INCLUDED

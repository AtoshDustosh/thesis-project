#include "SAM.h"



SAM* constructSAM(char* filePath, char* refName, uint64_t refLength) {
    SAM* sam = (SAM*)malloc(sizeof(SAM));
    if(sam == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }

    FILE* fpointer = fopen(filePath, "w");

    sam->filePath = filePath;
    if(fpointer == NULL){
        printf("ERROR: cannot create new file under file path - %s\n", filePath);
        exit(EXIT_FAILURE);
    }
    sam->fpointer = fpointer;
    sam->lastRead = NULL;
    sam->newRead = NULL;
    fprintf(sam->fpointer, "@HD VN:version\n");
    fprintf(sam->fpointer, "@SQ SN:%s LN:%"PRIu64"\n", refName, refLength);
    return sam;
}


void fprintReadUsingSAM(Read* read, SAM* sam) {
    if(sam->newRead == NULL){
        sam->newRead = read;
    } else if(sam->newRead != NULL && sam->lastRead == NULL) {
        sam->lastRead = sam->newRead;
        sam->newRead = read;
    } else {
        clearRead(sam->lastRead);
        sam->lastRead = sam->newRead;
        sam->newRead = read;
    }
    fprintf(sam->fpointer, "%s\t", sam->newRead->QNAME);
    fprintf(sam->fpointer, "%d\t", sam->newRead->FLAG);
    fprintf(sam->fpointer, "%s\t", sam->newRead->RNAME);
    fprintf(sam->fpointer, "%d\t", sam->newRead->POS);
    fprintf(sam->fpointer, "%f\t", sam->newRead->MAPQ);
    fprintf(sam->fpointer, "%s\t", sam->newRead->CIGAR);
    fprintf(sam->fpointer, "%s\t", sam->newRead->RNEXT);
    fprintf(sam->fpointer, "%d\t", sam->newRead->PNEXT);
    fprintf(sam->fpointer, "%d\t", sam->newRead->TLEN);
    fprintf(sam->fpointer, "%d\t", sam->newRead->FLAG);
    fprintf(sam->fpointer, "%s\t", sam->newRead->SEQ);
    fprintf(sam->fpointer, "%s\t", sam->newRead->QUAL);
    fprintf(sam->fpointer, "\n");

}

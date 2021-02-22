#include "Read.h"

#include <string.h>

void initRead(Read *read) {
    if(read == NULL) {
        printf("ERROR: null pointer occurred when initializing a read. \n");
        exit(EXIT_FAILURE);
    }
    read->QNAME = (char*)malloc(sizeof(char) * BUFSIZ);
    if(read->QNAME == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    read->QNAME[0] = '\0';
    read->FLAG = 0x0;
    read->RNAME = (char*)malloc(sizeof(char) * BUFSIZ);
    if(read->RNAME == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    read->RNAME[0] = '\0';
    read->POS = 0;
    read->MAPQ = 0;
    read->CIGAR = (char*)malloc(sizeof(char) * BUFSIZ);
    if(read->CIGAR == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    read->CIGAR[0] = '\0';
    read->RNEXT = (char*)malloc(sizeof(char) * BUFSIZ);
    if(read->RNEXT == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    read->RNEXT[0] = '\0';
    read->PNEXT = 0;
    read->TLEN = 0;
    read->SEQ = (char*)malloc(sizeof(char) * BUFSIZ);
    if(read->SEQ == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    read->SEQ[0] = '\0';
    read->QUAL = (char*)malloc(sizeof(char) * BUFSIZ);
    if(read->QUAL == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }
    read->QUAL[0] = '\0';
}




void printRead(Read* read) {
    if(read == NULL){
        printf("ERROR: null pointer occurs when printing a read.\n");
        exit(EXIT_FAILURE);
    }
    printf("QNAME\t%s\n", read->QNAME);
    printf("FLAG\t0x%x\n", read->FLAG);
    printf("RNAME\t%s\n", read->RNAME);
    printf("POS\t%d\n", read->POS);
    printf("MAPQ\t%f\n", read->MAPQ);
    printf("CIGAR\t%s\n", read->CIGAR);
    printf("RNEXT\t%s\n", read->RNEXT);
    printf("PNEXT\t%d\n", read->PNEXT);
    printf("TLEN\t%d\n", read->TLEN);
    printf("SEQ\t%s\n", read->SEQ);
    printf("QUAL\t%s\n", read->QUAL);
    printf("\n");

}

Read* copyRead(Read* read) {
    if(read == NULL){
        printf("ERROR: null pointer occurs when copying a read.\n");
        exit(EXIT_FAILURE);
    }
    Read* copiedRead = (Read*)malloc(sizeof(Read));
    if(copiedRead == NULL) {
        printf("ERROR: System memory not enough. \n");
        exit(EXIT_FAILURE);
    }

    initRead(copiedRead);
    strcpy(copiedRead->QNAME, read->QNAME);
    copiedRead->FLAG = read->FLAG;
    strcpy(copiedRead->RNAME, read->RNAME);
    copiedRead->POS = read->POS;
    copiedRead->MAPQ = read->MAPQ;
    strcpy(copiedRead->CIGAR, read->CIGAR);
    strcpy(copiedRead->RNEXT, read->RNEXT);
    copiedRead->PNEXT = read->PNEXT;
    copiedRead->TLEN = read->TLEN;
    strcpy(copiedRead->SEQ, read->SEQ);
    strcpy(copiedRead->QUAL, read->QUAL);
    return copiedRead;
}

void clearRead(Read* read) {
    if(read == NULL){
        printf("ERROR: null pointer occurs when clearing a read.\n");
        exit(EXIT_FAILURE);
    }
    free(read->QNAME);
    free(read->RNAME);
    free(read->CIGAR);
    free(read->RNEXT);
    free(read->SEQ);
    free(read->QUAL);
    free(read);
}









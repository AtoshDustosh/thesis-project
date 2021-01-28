#include "PairendAligner.h"






void loadSNAPStructure(SNAP* snap, PairendAligner* pairendAligner) {
    pairendAligner->snap = snap;
}


void initPairendAligner(char* fastqFilePath1, char* fastqFilePath2,
                        PairendAligner* pairendAligner) {
    pairendAligner->fastqFilePath1 = fastqFilePath1;
    pairendAligner->fastqFilePath2 = fastqFilePath2;

    /**
     * \todo integrate operations of FileOperation.h and main.c together.
     * They can be used here.
     *
     */
}

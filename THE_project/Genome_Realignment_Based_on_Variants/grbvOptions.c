#include "grbvOptions.h"


FileList* designatedFiles(Options *opts){
  int designatedFileCount = 0;
  if(getFaFile(opts) != NULL){
    designatedFileCount++;
  }
  if(getFastqFile(opts) != NULL){
    designatedFileCount++;
  }
  if(getSamFile(opts) != NULL){
    designatedFileCount++;
  }
  if(getvcfFile(opts) != NULL){
    designatedFileCount++;
  }

  if(designatedFileCount == 0){
    printf("Not designating any input file. Exiting the program ...\n");
    exit(EXIT_SUCCESS);
  }

  FileList* fileList = (FileList*)malloc(sizeof(FileList));
  fileList->paths = (char**)malloc(sizeof(char*) * designatedFileCount);
  fileList->count = designatedFileCount;
  int iteratorTemp = 0;
  if(getFaFile(opts) != NULL){
    fileList->paths[iteratorTemp] = opts->faFile;
    iteratorTemp++;
  }
  if(getFastqFile(opts) != NULL){
    fileList->paths[iteratorTemp] = opts->fastqFile;
    iteratorTemp++;
  }
  if(getSamFile(opts) != NULL){
    fileList->paths[iteratorTemp] = opts->samFile;
    iteratorTemp++;
  }
  if(getvcfFile(opts) != NULL){
    fileList->paths[iteratorTemp] = opts->vcfFile;
    iteratorTemp++;
  }

  return fileList;
}

void destroyFileList(FileList* fl){
  free(fl->paths);
  free(fl);
}
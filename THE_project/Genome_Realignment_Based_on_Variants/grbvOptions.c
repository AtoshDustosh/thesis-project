#include "grbvOptions.h"


FileList* designatedFiles(Options *opts){
  int designatedFileCount = 0;
  if(opts->faFile != NULL){
    designatedFileCount++;
  }
  if(opts->fastqFile != NULL){
    designatedFileCount++;
  }
  if(opts->samFile != NULL){
    designatedFileCount++;
  }
  if(opts->vcfFile != NULL){
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
  if(opts->faFile != NULL){
    fileList->paths[iteratorTemp] = opts->faFile;
    iteratorTemp++;
  }
  if(opts->fastqFile != NULL){
    fileList->paths[iteratorTemp] = opts->fastqFile;
    iteratorTemp++;
  }
  if(opts->samFile != NULL){
    fileList->paths[iteratorTemp] = opts->samFile;
    iteratorTemp++;
  }
  if(opts->vcfFile != NULL){
    fileList->paths[iteratorTemp] = opts->vcfFile;
    iteratorTemp++;
  }

  return fileList;
}

void destroyFileList(FileList* fl){
  free(fl->paths);
  free(fl);
}
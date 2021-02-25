#include <stdio.h>
#include <stdlib.h>


int main(int argc, char** argv){
	FILE* tmpFile = fopen("testFile.txt", "w");
   	
   	fprintf(tmpFile, "This is testing for fprintf...\n");
   	fputs("This is testing for fputs...\n", tmpFile);
	
	fclose(tmpFile);

	return 0;
}
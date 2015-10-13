#include "argon.h"
#include <stdio.h>
#include <stdlib.h>

int LoadData(char *input_file_name){

	int n;
	FILE *input_file;

	input_file=fopen(input_file_name,"r");
  
    if( input_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }
 
   fscanf(input_file, "%d,",&n);
   fclose(input_file);
   free(input_file);
   return n;
   
   }
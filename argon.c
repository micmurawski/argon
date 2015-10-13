#include "argon.h"
#include <stdio.h>
#include <stdlib.h>

int loadData(char *input_file_name){

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

 void saveData(char *output_file_name){

 }

 inline double *addArray(double array1[], double array2[], int n){

 	double *NewArray=(double*)malloc(n * sizeof(double)); 
 	for(int i=0;i<n;i++){
 		NewArray[i]=array1[i]+array2[i];
 	}

 	return NewArray;

 }


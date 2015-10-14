#include "argon.h"
#include <stdio.h>
#include <stdlib.h>

void loadData(char *input_file_name,int *data){

	FILE *input_file;

	input_file=fopen(input_file_name,"r");
  
    if( input_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }
 
   fscanf(input_file, "%d,",data);
   fclose(input_file);

   
   }

 void saveData(char *output_file_name,double data[][3], int N){

   FILE *output_file;

   output_file=fopen(output_file_name,"w");

    if( output_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }

   fprintf(output_file, "%d\n", N);
   for(int i=0;i<N;i++){
      fprintf(output_file, "atom%d %f %f %f\n",i+1,data[i][0],data[i][1],data[i][2]);
   }
   fclose(output_file);


 }

void addArray(double *result,double *array1, double *array2, int n){

   for(int i=0;i<3;i++){
      result[i]=array1[i]+array2[i];
   }

 }

 void multiplyArray(double *result,double *array1,double number,int n){

   for(int i=0;i<n;i++){
      result[i]=array1[i]*number;
   }

 }

void saveMomentum(char *output_file_name,double data[][3], int N){

   FILE *output_file;

   output_file=fopen(output_file_name,"w");

    if( output_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }

   for(int i=0;i<N;i++){
      fprintf(output_file, "%f \t %f \t %f\n",data[i][0],data[i][1],data[i][2]);
   }
   fclose(output_file);


}

double dotProduct(double *array1,double *array2, int n){
   double result;
   for(int i=0;i<n;i++){
      result+=array1[i]*array2[i];
   }
   return result;
}

double subtractArray(double *result,double *array1,double *array2, int n){
   
   for(int i=0;i<n;i++){
      result[i]=array1[i]-array2[i];
   }
}
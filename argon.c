#include "argon.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





void loadData(char *input_file_name,int intParam[],double doubleParam[]){

	FILE *input_file;

	input_file=fopen(input_file_name,"r");
  
    if( input_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }
 
   fscanf(input_file, "%d,",&intParam[0]);
   fscanf(input_file, "%lf,",&doubleParam[0]);
   fscanf(input_file, "%lf,",&doubleParam[1]);
   fscanf(input_file, "%lf,",&doubleParam[2]);
   fscanf(input_file, "%lf,",&doubleParam[3]);
   fscanf(input_file, "%lf,",&doubleParam[4]);
   fscanf(input_file, "%lf,",&doubleParam[5]);
   fscanf(input_file, "%lf,",&doubleParam[6]);
   fscanf(input_file, "%lf,",&doubleParam[7]);
   fscanf(input_file, "%d,",&intParam[1]);
   fscanf(input_file, "%d,",&intParam[2]);
   fscanf(input_file, "%d,",&intParam[3]);
   fscanf(input_file, "%d,",&intParam[4]);
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

double Norm(double *array1,int n){
   double result;
   for(int i=0;i<n;i++) result+=array1[i]*array1[i];
   return result;
}


void subtractArray(double *result,double *array1,double *array2, int n){
   
   for(int i=0;i<n;i++){
      result[i]=array1[i]-array2[i];
   }
}

double calculatePotential(double *array1, double *array2, int n, double eps, double R){

   double r;
   for(int i=0;i<n;i++){
      r+=pow(array1[i]-array2[i],2);
   }
   r=sqrt(r);

   return eps*(pow(R/r,12)-2*pow(R/r,6));

}

void calculateForce(double *result, double *array1, double *array2, int n, double eps, double R){

   double r;
   for(int i=0;i<n;i++){
      r+=pow(array1[i]-array2[i],2);
   }
   r=sqrt(r);

   subtractArray(result,array1,array2,n);

   multiplyArray(result,result,
      12*eps*(pow(R/r,12)-pow(R/r,6))/(r*r),
      3);

}
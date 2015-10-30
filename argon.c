#include "argon.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>





inline void loadData(char *input_file_name,int intParam[],double doubleParam[]){

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

 inline void saveData(char *output_file_name,double data[][3], int N){

   FILE *output_file;

   output_file=fopen(output_file_name,"w");

    if( output_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }

   fprintf(output_file, "%d\n", N);
   fprintf(output_file, "comment\n");
   for(int i=0;i<N;i++){
      fprintf(output_file, "atom%d\t%f\t%f\t%f\n",i,data[i][0],data[i][1],data[i][2]);
   }
   fclose(output_file);


 }


inline void savePositions(FILE *output_file,double data[][3],int N){

   for(int i=0;i<N;i++){
      fprintf(output_file, "atom%d\t%f\t%f\t%f\n",i,data[i][0],data[i][1],data[i][2]);
   }
   fprintf(output_file, "\n");

}


inline void saveMomentum(char *output_file_name,double data[][3], int N){

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


inline double distance(double *array1,double *array2, int n){
   double result;
   for(int i;i<n;i++) result+=(array1[i]-array2[i])*(array1[i]-array2[i]);
   return sqrt(result);
}

inline double temperature(double data[][3], double mass,double k, int n){
   double T=0;
   for(int i=0;i<n;i++) T+=((data[i][0]*data[i][0])+(data[i][1]*data[i][1])+(data[i][2]*data[i][2]));
   T/=(3.0*k*n*mass);
   return T;
}

inline double pressure(double data[][3], double L, int N){
   double result;
   for(int i;i<N;i++)result+=sqrt((data[i][0]*data[i][0])+(data[i][1]*data[i][1])+(data[i][2]*data[i][2]));
   return result/(4*3.1415926535*L*L);
}

inline double sum(double data[], int N){
   double result;
   for (int i;i<N;i++) result+=data[i];
   return result;
}

inline double kineticEnergy(double data[][3],double mass, int N){
   double result=0;
   for(int i=0;i<N;i++)result+=((data[i][0]*data[i][0])+(data[i][1]*data[i][1])+(data[i][2]*data[i][2]));
   return result/(2*mass);
}
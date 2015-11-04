#include "argon.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline void loadData(char *input_file_name,struct Parameters *p){

   FILE *input_file;

   input_file=fopen(input_file_name,"r");
  
    if( input_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }
 
   fscanf(input_file, "%d,",&p->n);
   fscanf(input_file, "%lf,",&p->m);
   fscanf(input_file, "%lf,",&p->eps);
   fscanf(input_file, "%lf,",&p->R);
   fscanf(input_file, "%lf,",&p->f);
   fscanf(input_file, "%lf,",&p->L);
   fscanf(input_file, "%lf,",&p->A);
   fscanf(input_file, "%lf,",&p->T);
   fscanf(input_file, "%lf,",&p->TAU);
   fscanf(input_file, "%d,",&p->s_0);
   fscanf(input_file, "%d,",&p->s_d);
   fscanf(input_file, "%d,",&p->s_out);
   fscanf(input_file, "%d,",&p->s_xyz);
   fclose(input_file);
   
   }


inline double temperature(double *px,double *py,double *pz, double mass,double k, int n){
   double T=0;
   for(int i=0;i<n;i++) T+=((px[i]*px[i])+(py[i]*py[i])+(pz[i]*pz[i]));
   T/=(3.0*k*n*mass);
   return T;
}

inline double pressure(double *px,double *py,double *pz, double L, int N){
   double result=0;
   for(int i=0;i<N;i++)result+=sqrt((px[i]*px[i])+(py[i]*py[i])+(pz[i]*pz[i]));
   return result/(4*3.1415926535*L*L);
}

inline double kineticEnergy(double *px,double *py,double *pz,double mass, int N){
   double result=0;
   for(int i=0;i<N;i++)result+=((px[i]*px[i])+(py[i]*py[i])+(pz[i]*pz[i]));
   return result/(2*mass);
}

inline double sumArray(double *v, int N){
   double result=0;
   for(int i=0;i<N;i++)result+=v[i];
   return result;
}

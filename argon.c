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
   fscanf(input_file, "%Lf,",&p->m);
   fscanf(input_file, "%Lf,",&p->eps);
   fscanf(input_file, "%Lf,",&p->R);
   fscanf(input_file, "%Lf,",&p->f);
   fscanf(input_file, "%Lf,",&p->L);
   fscanf(input_file, "%Lf,",&p->A);
   fscanf(input_file, "%Lf,",&p->T);
   fscanf(input_file, "%Lf,",&p->TAU);
   fscanf(input_file, "%d,",&p->s_0);
   fscanf(input_file, "%d,",&p->s_d);
   fscanf(input_file, "%d,",&p->s_out);
   fscanf(input_file, "%d,",&p->s_xyz);
   fclose(input_file);
   
   }


inline long double temperature(long double *px,long double *py,long double *pz, long double mass,long double k, int n){
   long double T=0;
   for(int i=0;i<n;i++) T+=((px[i]*px[i])+(py[i]*py[i])+(pz[i]*pz[i]));
   T/=(3.0*k*n*mass);
   return T;
}

inline long double pressure(long double *px,long double *py,long double *pz, long double L, int N){
   long double result=0;
   for(int i=0;i<N;i++)result+=sqrtl((px[i]*px[i])+(py[i]*py[i])+(pz[i]*pz[i]));
   return result/(4*3.1415926535*L*L);
}

inline long double kineticEnergy(long double *px,long double *py,long double *pz,long double mass, int N){
   long double result=0;
   for(int i=0;i<N;i++)result+=((px[i]*px[i])+(py[i]*py[i])+(pz[i]*pz[i]));
   return result/(2*mass);
}

inline long double sumArray(long double *v, int N){
   long double result=0;
   for(int i=0;i<N;i++)result+=v[i];
   return result;
}

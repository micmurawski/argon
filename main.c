#include <stdio.h>
#include <stdlib.h>
#include "argon.h"
#include <math.h>
#include <time.h>
#include <string.h>

#define INPUT  argv[1]
#define OUTPUT argv[2]
#define OUTPUT_xyz argv[3]

#define PI 3.1415926535
#define K   8.31e-3

#define n parameters.n
#define M parameters.m
#define EPS parameters.eps
#define R parameters.R
#define f parameters.f
#define L parameters.L
#define A parameters.A
#define T0 parameters.T
#define TAU parameters.TAU
#define s_0 parameters.s_0
#define s_d parameters.s_d
#define s_out parameters.s_out
#define s_xyz parameters.s_xyz

int main(int argc, char *argv[]){

//if(strcmp(argv[2],"-")!=0){
//  FILE *output_file;
//}
FILE *output_xyz_file;

struct Parameters parameters;
loadData(INPUT,&parameters);
double apx,apy,apz;
double sumT=0,sumK=0,sumV=0,sumP=0,Vt,Pt,Kt,Tt;
double random_number;
int i0,i1,i2,i;
const int N=n*n*n;

double* px = (double *)malloc(N*sizeof(double));
double* py = (double *)malloc(N*sizeof(double));
double* pz = (double *)malloc(N*sizeof(double));

double* Fx = (double *)calloc(N,sizeof(double));
double* Fy = (double *)calloc(N,sizeof(double));
double* Fz = (double *)calloc(N,sizeof(double));

double* x = (double *)calloc(N,sizeof(double));
double* y = (double *)calloc(N,sizeof(double));
double* z = (double *)calloc(N,sizeof(double));

double* V = (double *)calloc(N,sizeof(double));

srand(time(NULL));
   
   printf("INPUT FILE: %s\n", INPUT);
   printf("OUTPUT FILE:%s\n", OUTPUT_xyz);
   printf("ATOM NUMBER:%d\n",N);
   printf("ATOM MASS:  %lf u\n",M);
   printf("EPS:        %lf KJ/mol\n",EPS);
   printf("R:          %lf nm\n",R);
   printf("f:          %lf\n",f);
   printf("L:          %lf nm\n",L);
   printf("A:          %lf nm\n",A);
   printf("T:          %lf K\n",T0);
   printf("TAU:        %lf ps\n",TAU);
   printf("S_0:        %d\n",s_0);
   printf("S_D:        %d\n",s_d);
   printf("S_OUT:      %d\n",s_out);
   printf("S_XYZ:      %d\n",s_xyz);
   

//generating initial positions of atoms
   for(i0=0;i0<n;i0++){
      for(i1=0;i1<n;i1++){
         for(i2=0;i2<n;i2++){

            i=i0+i1*n+i2*n*n;

            x[i]+=(double)(i0-(n-1)/2)*A;

            x[i]+=(double)(i1-(n-1)/2)*A/2;
            y[i]+=(double)(i1-(n-1)/2)*A*sqrt(3)/2;
            
            x[i]+=(double)(i2-(n-1)/2)*A/2;
            y[i]+=(double)(i2-(n-1)/2)*A*sqrt(3)/6;
            z[i]+=(double)(i2-(n-1)/2)*A*sqrt(2.0/3.0);

         }
         
      }
            }

   output_xyz_file=fopen(OUTPUT_xyz,"w");

    if( output_xyz_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }

   for(int ii=0;ii<N;ii++){
     
     random_number=(((double)rand()+1)/((double)RAND_MAX+1)); //  random number (0;1]
     px[ii]=sqrt(-2*K*T0*M*logl(random_number)); // generating initial momentum form maxwell-boltzman distribution
     random_number=((double)rand()/RAND_MAX); //  random number [0;1]
     px[ii]*=cos(2*PI*random_number); // generating sign of momentum according to box-muller transformation

     random_number=(((double)rand()+1)/((double)RAND_MAX+1));
     py[ii]=sqrt(-2*K*T0*M*logl(random_number));
     random_number=((double)rand()/RAND_MAX);
     py[ii]*=cos(2*PI*random_number); 

     random_number=(((double)rand()+1)/((double)RAND_MAX+1)); 
     pz[ii]=sqrt(-2*K*T0*M*logl(random_number)); 
     random_number=((double)rand()/RAND_MAX);
     pz[ii]*=cos(2*PI*random_number); 


   }

    apx=sumArray(px,N)/N;
    apy=sumArray(py,N)/N;
    apz=sumArray(pz,N)/N;

   for(int ii=0;ii<N;ii++){
    px[ii]-=apx;
    py[ii]-=apy;
    pz[ii]-=apz;
   }


  algorytm2(x,y,z,Fx,Fy,Fz,V,EPS,R,f,L,N);

   printf("AVERAGE MOMENTUM %lf %lf %lf\n",sumArray(px,N)/N,sumArray(py,N)/N,sumArray(pz,N)/N);
   printf("CALCULATED TEMPERATURE %lf \n", temperature(px,py,pz,M,K,N) );
   printf("CALCULATED KIN ENERGY %lf \n", kineticEnergy(px,py,pz,M,N));
   printf("CALCULATED POT ENERGY %lf \n", sumArray(V,N));
   printf("CALCULATED TOTAL ENERGY %lf \n", kineticEnergy(px,py,pz,M,N)+sumArray(V,N));
   printf("CALCULATED PRESSURE %lf \n", pressure(px,py,pz,L,N));
   getchar();
   printf("TIME [ps]\tE_TOT [kJ/mol]\tE_POT [kJ/mol]\tTEMP. [K]\tPRES. [kJ/mol/nm^3]\n");

   //ILOSC KROKOW
   for(int ss=0;ss<s_d+s_0+1;ss++){

  for(int ii=0;ii<N;ii++){

  px[ii]+=0.5*Fx[ii]*TAU;
  py[ii]+=0.5*Fy[ii]*TAU;
  pz[ii]+=0.5*Fz[ii]*TAU;

  x[ii]+=px[ii]*TAU/M;
  y[ii]+=py[ii]*TAU/M;
  z[ii]+=pz[ii]*TAU/M;
}

  algorytm2(x,y,z,Fx,Fy,Fz,V,EPS,R,f,L,N);

  Vt=sumArray(V,N);
  if(ss>s_0) sumV+=Vt;
  Pt=pressure(px,py,pz,L,N);
  sumP+=Pt;

   for(int ii=0;ii<N;ii++){

   px[ii]+=0.5*Fx[ii]*TAU;
   py[ii]+=0.5*Fy[ii]*TAU;
   pz[ii]+=0.5*Fz[ii]*TAU;
}

Kt=kineticEnergy(px,py,pz,M,N);
if(ss>s_0) sumK+=Kt;
Tt=temperature(px,py,pz,M,K,N);
if(ss>s_0) sumT+=Tt;

//OBLICZENIE NOWYCH PEDOW I NOWYCH POLOZEN
if(ss==s_0 && s_0 >0){
printf("KONIEC TERMALIZACJI\n");
//sumV=0;
//sumP=0;
//sumK=0;
//sumT=0;
}

if(ss>s_0){

if(ss%s_xyz==0){
   fprintf(output_xyz_file, "%d\n",N);
   fprintf(output_xyz_file, "ss=%lf\n",(ss-s_0)*TAU);
  for(int i=0; i<N;i++)fprintf(output_xyz_file, "atom%d\t%lf\t%lf\t%lf\n",i,x[i],y[i],z[i]);
}

if(ss%s_out==0){
  printf("%lf\t%lf\t%lf\t%lf\t%lf\t\n",(ss-s_0)*TAU,Kt+Vt,Vt,Tt,Pt);
}
  

}
}
    printf("AVERAGE\t%lf\t%lf\t%lf\t%lf\t\n",(sumK+sumV)/s_d,sumV/s_d,sumT/s_d,sumP/s_d);
    free(Fx);
    free(Fy);
    free(Fz);

    free(px);
    free(py);
    free(pz);

    free(x);
    free(y);
    free(z);

    fclose(output_xyz_file);

    return 0;

}
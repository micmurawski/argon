#include <stdio.h>
#include <stdlib.h>
#include "argon.h"
#include <math.h>
#include <time.h>
#include <string.h>

//#define Parameters struct Parameters
#define INPUT  argv[1]
#define OUTPUT argv[2]
#define OUTPUT_xyz argv[3]

#define PI 3.1415926535
#define K   8.31e-5

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

struct Parameters parameters;
loadData(INPUT,&parameters);
FILE *output_xyz_file;

double rji,rj;
int i0,i1,i2,i;
double random_number;
const int N=n*n*n;

double* px = (double *)malloc(N*sizeof(double));
double* py = (double *)malloc(N*sizeof(double));
double* pz = (double *)malloc(N*sizeof(double));

double* Fx = (double *)calloc(N,sizeof(double));
double* Fy = (double *)calloc(N,sizeof(double));
double* Fz = (double *)calloc(N,sizeof(double));

double* rx = (double *)calloc(N,sizeof(double));
double* ry = (double *)calloc(N,sizeof(double));
double* rz = (double *)calloc(N,sizeof(double));

double* V = (double *)calloc(N,sizeof(double));

srand(time(NULL));
   

   printf("INPUT FILE: %s\n", INPUT);
   printf("OUTPUT FILE: %s\n", OUTPUT_xyz);
   printf("ATOM NUMBER: %d\n",N);
   printf("TEMPERATURE: %f\n",T0);
   printf("A:           %f\n",A);
   printf("ATOM MASS:   %f\n",M);
   printf("EPS:         %f\n",EPS);
   printf("R:           %f\n",R);
   printf("TAU:           %f\n",TAU);
   printf("f:           %f\n",f);
   printf("L:           %f\n",L);

//generating initial positions of atoms
   for(i0=0;i0<n;i0++){
      for(i1=0;i1<n;i1++){
         for(i2=0;i2<n;i2++){

            i=i0+i1*n+i2*n*n;

            rx[i]+=(double)(i0-(n-1)/2)*A;

            rx[i]+=(double)(i1-(n-1)/2)*A/2;
            ry[i]+=(double)(i1-(n-1)/2)*A*sqrt(3)/2;
            
            rx[i]+=(double)(i2-(n-1)/2)*A/2;
            ry[i]+=(double)(i2-(n-1)/2)*A*sqrt(3)/6;
            rz[i]+=(double)(i2-(n-1)/2)*A*sqrt(2.0/3.0);

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
     px[ii]=sqrt(-2*K*T0*M*log(random_number)); // generating initial momentum form maxwell-boltzman distribution
     random_number=((double)rand()/RAND_MAX); //  random number [0;1]
     px[ii]*=cos(2*PI*random_number); // generating sign of momentum according to box-muller transformation

     random_number=(((double)rand()+1)/((double)RAND_MAX+1));
     py[ii]=sqrt(-2*K*T0*M*log(random_number));
     random_number=((double)rand()/RAND_MAX);
     py[ii]*=cos(2*PI*random_number); 

     random_number=(((double)rand()+1)/((double)RAND_MAX+1)); 
     pz[ii]=sqrt(-2*K*T0*M*log(random_number)); 
     random_number=((double)rand()/RAND_MAX);
     pz[ii]*=cos(2*PI*random_number); 


   }

   double apx=sumArray(px,N)/N;
   double apy=sumArray(py,N)/N;
   double apz=sumArray(pz,N)/N;
   for(int ii=0;ii<N;ii++){
    px[ii]-=apx;
    py[ii]-=apy;
    pz[ii]-=apz;
   }

   printf("AVERAGE MOMENTUM %f %f %f\n",sumArray(px,N)/N,sumArray(py,N)/N,sumArray(pz,N)/N);
   printf("CALCULATED TEMPERATURE %f \n", temperature(px,py,pz,M,K,N) );
   printf("CALCULATED KIN ENERGY %f \n", kineticEnergy(px,py,pz,M,N) );
   getchar();

   
   fprintf(output_xyz_file, "%d\n",N);
   fprintf(output_xyz_file, "ss=%d\n",0);
   for(int i=0; i<N;i++)fprintf(output_xyz_file, "atom%d\t%f\t%f\t%f\n",i,rx[i],ry[i],rz[i]);

    for(int jj=0;jj<N;jj++){
   for(int ii=jj;ii<N;ii++){
   
   // Rji=|Rj-Ri|
    if(ii!=jj){

   rji=(rx[jj]-rx[ii])*(rx[jj]-rx[ii]);
   rji+=(ry[jj]-ry[ii])*(ry[jj]-ry[ii]);
   rji+=(rz[jj]-rz[ii])*(rz[jj]-rz[ii]);
   rji=sqrt(rji);


   V[jj]+=EPS*(pow(R/rji,12)-2*pow(R/rji,6));
   V[ii]+=EPS*(pow(R/rji,12)-2*pow(R/rji,6));
   


   Fx[jj]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rx[jj]-rx[ii])/(rji*rji);

   Fy[jj]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(ry[jj]-ry[ii])/(rji*rji);

   Fz[jj]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rz[jj]-rz[ii])/(rji*rji);

   Fx[ii]-=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rx[jj]-rx[ii])/(rji*rji);

   Fy[ii]-=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(ry[jj]-ry[ii])/(rji*rji);

   Fz[ii]-=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rz[jj]-rz[ii])/(rji*rji);


}
 }

   rj=rx[jj]*rx[jj];
   rj+=ry[jj]*ry[jj];
   rj+=rz[jj]*rz[jj];
   rj=sqrt(rj);

   if(rj>=L){

  V[jj]+=0.5*f*(rj-L)*(rj-L);

   Fx[jj]+=f*(L-rj)*rx[jj]/rj;

   Fy[jj]+=f*(L-rj)*ry[jj]/rj;

   Fz[jj]+=f*(L-rj)*rz[jj]/rj;

   }

   }


   //ILOSC KROKOW
   for(int ss=0;ss<2000;ss++){

  for(int ii=0;ii<N;ii++){

   px[ii]+=0.5*Fx[ii]*TAU;
   py[ii]+=0.5*Fy[ii]*TAU;
   pz[ii]+=0.5*Fz[ii]*TAU;

   rx[ii]+=px[ii]*TAU/M;
   ry[ii]+=py[ii]*TAU/M;
   rz[ii]+=pz[ii]*TAU/M;

}



   for(int jj=0;jj<N;jj++){
   for(int ii=jj;ii<N;ii++){
   
   // Rji=|Rj-Ri|
    if(ii!=jj){

   rji=(rx[jj]-rx[ii])*(rx[jj]-rx[ii]);
   rji+=(ry[jj]-ry[ii])*(ry[jj]-ry[ii]);
   rji+=(rz[jj]-rz[ii])*(rz[jj]-rz[ii]);
   rji=sqrt(rji);


   V[jj]+=EPS*(pow(R/rji,12)-2*pow(R/rji,6));
   V[ii]+=EPS*(pow(R/rji,12)-2*pow(R/rji,6));
   


   Fx[jj]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rx[jj]-rx[ii])/(rji*rji);

   Fy[jj]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(ry[jj]-ry[ii])/(rji*rji);

   Fz[jj]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rz[jj]-rz[ii])/(rji*rji);

   Fx[ii]-=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rx[jj]-rx[ii])/(rji*rji);

   Fy[ii]-=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(ry[jj]-ry[ii])/(rji*rji);

   Fz[ii]-=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(rz[jj]-rz[ii])/(rji*rji);


}
 }

   rj=rx[jj]*rx[jj];
   rj+=ry[jj]*ry[jj];
   rj+=rz[jj]*rz[jj];
   rj=sqrt(rj);

   if(rj>=L){

  V[jj]+=0.5*f*(rj-L)*(rj-L);

   Fx[jj]+=f*(L-rj)*rx[jj]/rj;

   Fy[jj]+=f*(L-rj)*ry[jj]/rj;

   Fz[jj]+=f*(L-rj)*rz[jj]/rj;

   }

   }

   for(int ii=0;ii<N;ii++){

   px[ii]+=0.5*Fx[ii]*TAU;
   py[ii]+=0.5*Fy[ii]*TAU;
   pz[ii]+=0.5*Fz[ii]*TAU;
}

//OBLICZENIE NOWYCH PEDOW I NOWYCH POLOZEN


//for(int k=0; k)
  memset(Fy,0,sizeof(double)*N);
  memset(Fx,0,sizeof(double)*N);
  memset(Fz,0,sizeof(double)*N);
  

  //printf("kasujemy!\n");
  //for(int ii=0;ii<N;ii++)printf("%f %f %f \n",Fx[ii],Fy[ii],Fz[ii]);

   fprintf(output_xyz_file, "%d\n",N);
   fprintf(output_xyz_file, "ss=%d\n",ss);
//ZAPISANIE POZYCJI W FUNKCJI KROKU CZASOWEGO

  for(int i=0; i<N;i++)fprintf(output_xyz_file, "atom%d\t%f\t%f\t%f\n",i,rx[i],ry[i],rz[i]);
  //for(int i=0;i<N;i++)printf("%f \n",V[i]);

  printf("%d ",ss);
  printf("%f ",temperature(px,py,pz,M,K,N));
  printf("%f ",sumArray(V,N));
  printf("%f \n",kineticEnergy(px,py,pz,M,N));
  memset(V,0,sizeof(double)*N);
   //for(int i=0;i<N;i++)printf("%f\n",V[i]);

  


   

}
    
    free(Fx);
    free(Fy);
    free(Fz);

    free(px);
    free(py);
    free(pz);

    free(rx);
    free(ry);
    free(rz);

    fclose(output_xyz_file);

    return 0;

}
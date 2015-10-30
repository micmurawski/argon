#include <stdio.h>
#include <stdlib.h>
#include "argon.h"
#include <math.h>
#include <time.h>

#define Parameters struct Parameters
#define INPUT  argv[1]
#define OUTPUT argv[2]
#define PI 3.1415926535
#define K   8.31e-5

#define n intParam[0]
#define M doubleParam[0]
#define EPS doubleParam[1]
#define R doubleParam[2]
#define f doubleParam[3]
#define L doubleParam[4]
#define A doubleParam[5]
#define T0 doubleParam[6]
#define TAU doubleParam[7]



int main(int argc, char *argv[]){

int intParam[5];
double doubleParam[8];
loadData(INPUT,intParam,doubleParam);
FILE *output_file;

double rji;
int i0,i1,i2,i;
double random_number;
double averageP[3]={0,0,0};


srand(time(NULL));
   
   printf("INPUT FILE: %s\n", INPUT);
   printf("OUTPUT FILE: %s\n", OUTPUT);
   const int N=n*n*n;
   double r[N][3],p[N][3],V[N], F[N][3];
   printf("ATOM NUMBER: %d\n",N);
   printf("TEMPERATURE: %f\n",T0);
   printf("AVERAGE ENERGY:%f\n",K*T0*N*1/2);
   printf("A:           %f\n",A);
   printf("ATOM MASS:   %f\n",M);
   printf("EPS:         %f\n",EPS);
   printf("R:           %f\n",R);
   printf("TAU:           %f\n",TAU);
   printf("f:           %f\n",f);
   printf("L:           %f\n",L);
   printf("B0=(%f, %f, %f)\n",A,0.0,0.0);
   printf("B1=(%f, %f, %f)\n",A/2.0,A*sqrt(3)/2,0.0);
   printf("B2=(%f, %f, %f)\n",A/2.0,A*sqrt(3)/6,A*sqrt(2.0/3.0));

   

//generating initial positions of atoms
   for(i0=0;i0<n;i0++){
      for(i1=0;i1<n;i1++){
         for(i2=0;i2<n;i2++){

            i=i0+i1*n+i2*n*n;

            r[i][0]+=(double)(i0-(n-1)/2)*A;

            r[i][0]+=(double)(i1-(n-1)/2)*A/2;
            r[i][1]+=(double)(i1-(n-1)/2)*A*sqrt(3)/2;
            
            r[i][0]+=(double)(i2-(n-1)/2)*A/2;
            r[i][1]+=(double)(i2-(n-1)/2)*A*sqrt(3)/6;
            r[i][2]+=(double)(i2-(n-1)/2)*A*sqrt(2.0/3.0);

         }
         
      }
            }
            



   output_file=fopen(OUTPUT,"w");

    if( output_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }
   fprintf(output_file, "%d\n",N);
   fprintf(output_file, "comment\n");

            //for(int ii=0;ii<n*n*n;ii++) printf("atom%d (%f, %f, %f)\n",ii,r[ii][0],r[ii][1],r[ii][2]);
            savePositions(output_file,r,N);

            

   for(int ii=0;ii<N;ii++){
     
     random_number=(((double)rand()+1)/((double)RAND_MAX+1)); //  random number (0;1]
     p[ii][0]=sqrt(-2*K*T0*M*log(random_number)); // generating initial momentum form maxwell-boltzman distribution
     random_number=((double)rand()/RAND_MAX); //  random number [0;1]
     p[ii][0]*=cos(2*PI*random_number); // generating sign of momentum according to box-muller transformation

     random_number=(((double)rand()+1)/((double)RAND_MAX+1));
     p[ii][1]=sqrt(-2*K*T0*M*log(random_number));
     random_number=((double)rand()/RAND_MAX);
     p[ii][1]*=cos(2*PI*random_number); 

     random_number=(((double)rand()+1)/((double)RAND_MAX+1)); 
     p[ii][2]=sqrt(-2*K*T0*M*log(random_number)); 
     random_number=((double)rand()/RAND_MAX);
     p[ii][2]*=cos(2*PI*random_number); 


   }

   printf("KINTEIC ENERGY: %f\n",kineticEnergy(p,M,N));

   saveMomentum("output2",p,N);
   printf("averageP= %f %f %f\n", sum(p[0],N),sum(p[1],N),sum(p[2],N));
   printf("CALCULATED TEMPERATURE %f \n", temperature(p,M,K,N) );
   getchar();

   //ILOSC KROKOW
   for(int ss=0;ss<1000;ss++){
   //printf("%f \n",temperature(p,M,K,N));

   //OBLICZENIE POTENCJALOW I SIL
   //saveMomentum("output3",F,N);

   memset(&F,0,sizeof(F));
   for(int jj;jj<N;jj++){
   for(int ii=jj;ii<N;ii++){
   

   // Rji=|Rj-Ri|
      if(ii!=jj){

   rji=distance(r[jj],r[ii],3);

   F[jj][0]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(r[jj][0]-r[ii][0])/(rji*rji);
   F[jj][1]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(r[jj][1]-r[ii][1])/(rji*rji);
   F[jj][2]+=12*EPS*(pow(R/rji,12)-pow(R/rji,6))*(r[jj][2]-r[ii][2])/(rji*rji);

   F[jj+ii][0]=-1*F[jj][0];
   F[jj+ii][1]=-1*F[jj][1];
   F[jj+ii][2]=-1*F[jj][2];
}
      }


   }
   
//OBLICZENIE NOWYCH PEDOW I NOWYCH POLOZEN
for(int ii=0;ii<N;ii++){
   
   p[ii][0]+=0.5*F[ii][0]*TAU;
   p[ii][1]+=0.5*F[ii][1]*TAU;
   p[ii][2]+=0.5*F[ii][2]*TAU;

   r[ii][0]+=p[ii][0]*TAU;
   r[ii][1]+=p[ii][1]*TAU;
   r[ii][2]+=p[ii][2]*TAU;

}

 fprintf(output_file, "%d\n",N);
   fprintf(output_file, "comment\n");
//ZAPISANIE POZYCJI W FUNKCJI KROKU CZASOWEGO
savePositions(output_file,r,N);
printf("%d %.1f %f\n",ss,temperature(p,M,K,N),kineticEnergy(p,M,N));

}
fclose(output_file);


}
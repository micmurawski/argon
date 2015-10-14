#include <stdio.h>
#include <stdlib.h>
#include "argon.h"
#include <math.h>
#include <time.h>
#define INPUT	argv[1]
#define OUTPUT	argv[2]
#define PI 3.1415926535



int main(int argc, char *argv[]){

const double A=1;
const double M=1;
const double T=100;
const double K=8.6173303e-5;
const double EPS=1;
const double R=1;

int n,i0,i1,i2,i;
float random_number;
double averageP[3]={0,0,0};
double v0[3]={0,0,0};
const double b0[3]= {A,0,0};
const double b1[3]= {A/2,A*sqrt(3)/2,0};
const double b2[3]= {A/2,A*sqrt(3)/6,A*sqrt(2.0/3.0)};

srand(time(NULL));
	
	printf("INPUT FILE: %s\n", INPUT);
	printf("OUTPUT FILE: %s\n", OUTPUT);
	loadData(INPUT,&n);

	int N=n*n*n;
	double r[N][3],p[N][3], V[N][3];
   printf("ATOM NUMBER: %d\n",N);
	printf("B0=(%f, %f, %f)\n",b0[0],b0[1],b0[2]);
	printf("B1=(%f, %f, %f)\n",b1[0],b1[1],b1[2]);
	printf("B2=(%f, %f, %f)\n",b2[0],b2[1],b2[2]);

   for(i0=0;i0<n;i0++){
   	for(i1=0;i1<n;i1++){
   		for(i2=0;i2<n;i2++){

   			i=i0+i1*n+i2*n*n;

   			memcpy(&v0, &b0, sizeof v0);
   			multiplyArray(v0,v0,(double)(i0-(n-1)/2),3);
   			addArray(r[i],r[i],v0,3);

   			memcpy(&v0, &b1, sizeof v0);
   			multiplyArray(v0,v0,(double)(i1-(n-1)/2),3);
   			addArray(r[i],r[i],v0,3);

   			memcpy(&v0, &b2, sizeof v0);
   			multiplyArray(v0,v0,(double)(i2-(n-1)/2),3);
   			addArray(r[i],r[i],v0,3);
   		}
   		
   	}
   			}
				
				for(int ii=0;ii<n*n*n;ii++) printf("atom%d (%f, %f, %f)\n",ii,r[ii][0],r[ii][1],r[ii][2]);
				saveData(OUTPUT,r,N);

   for(int ii=0;ii<N;ii++){
     
     random_number=(((double)rand()+1)/((double)RAND_MAX+1)); //  random number (0;1]

      p[ii][0]=sqrt(-2*K*T*M*log(random_number));
      random_number=((double)rand()/RAND_MAX);
      p[ii][0]*=cos(2*PI*random_number);
      //if(random_number>=0.5) p[ii][0]*=-1;
      //p[ii][0]*=cos(2*3.14*random_number);

      
      for(int ii=0;ii<N;ii++){
         addArray(averageP,averageP,p[ii],3);
         multiplyArray(averageP,averageP,((double)1/N),3);
      }

   }

   saveMomentum("output2",p,N);
   printf("averageP= %f %f %f\n", averageP[0],averageP[1],averageP[2]);




}
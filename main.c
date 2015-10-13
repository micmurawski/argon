#include <stdio.h>
#include <stdlib.h>
#include "argon.h"
#include <math.h>
const double a=1;

int i;

char *INPUT,*OUTPUT;
int main(int argc, char *argv[]){

int n,i0,i1,i2;
double v0[3]={0,0,0};
double b0[3]= {a,0,0};
double b1[3]= {a/2,a*sqrt(3)/2,0};
double b2[3]= {a/2,a*sqrt(3)/6,a*sqrt(2.0/3.0)};
	
	INPUT=argv[1];
	OUTPUT=argv[2];
	printf("INPUT FILE: %s\n", INPUT);
	printf("OUTPUT FILE: %s\n", OUTPUT);
	loadData(argv[1],&n);
	int N=n*n*n;
	double r[N][3],p[N][3];
   printf("liczba atomow: %d\n",N);
	printf("b0=(%f, %f, %f)\n",b0[0],b0[1],b0[2]);
	printf("b1=(%f, %f, %f)\n",b1[0],b1[1],b1[2]);
	printf("b2=(%f, %f, %f)\n",b2[0],b2[1],b2[2]);

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
   			//memcpy(&v0, &b0, sizeof v0);
   			//multiplyArray(v0,v0,(i0-(n-1)/2),3);
   			//printf("%f,%f,%f\n",v0[0],v0[1],v0[2]);
   			//printf("hej %d %d %d\n",i0,i1,i2);
   			//addArray(r[i],r[i],v0,3);

				
				for(int ii=0;ii<n*n*n;ii++) printf("atom%d  (%f, %f, %f)\n",ii,r[ii][0],r[ii][1],r[ii][2]);
				saveData(OUTPUT,r,N);




}
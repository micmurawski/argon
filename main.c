#include <stdio.h>
#include <stdlib.h>
#include "argon.h"
#include <math.h>
#include <time.h>

#define Parameters struct Parameters
#define INPUT	argv[1]
#define OUTPUT	argv[2]
#define PI 3.1415926535
#define K   8.31e-5
#define T0 doubleParam[6]
#define n intParam[0]
#define A doubleParam[5]
#define M doubleParam[0]
#define EPS doubleParam[1]
#define R doubleParam[2]
#define TAU doubleParam[7]
#define F doubleParam[3]
#define L doubleParam[2]



int main(int argc, char *argv[]){

int intParam[5];
double doubleParam[8];
loadData(INPUT,intParam,doubleParam);
FILE *output_file;

double b0[3]= {A,0,0};
double b1[3]= {A/2,A*sqrt(3)/2,0};
double b2[3]= {A/2,A*sqrt(3)/6,A*sqrt(2.0/3.0)};
double force[3];
double sumforce[3];

int i0,i1,i2,i;
double random_number;
double averageP[3]={0,0,0};
double v0[3]={0,0,0};


srand(time(NULL));
	
	printf("INPUT FILE: %s\n", INPUT);
	printf("OUTPUT FILE: %s\n", OUTPUT);
	const int N=n*n*n;
	double r[N][3],p[N][3],V[N];
   printf("ATOM NUMBER: %d\n",N);
   printf("TEMPERATURE: %f\n",T0);
   printf("AVERAGE ENERGY:%f\n",K*T0*N*1/2);
   printf("A:           %f\n",A);
   printf("ATOM MASS:   %f\n",M);
   printf("EPS:         %f\n",EPS);
   printf("R:           %f\n",R);
   printf("TAU:           %f\n",TAU);
	printf("B0=(%f, %f, %f)\n",b0[0],b0[1],b0[2]);
	printf("B1=(%f, %f, %f)\n",b1[0],b1[1],b1[2]);
	printf("B2=(%f, %f, %f)\n",b2[0],b2[1],b2[2]);

//generating initial positions of atoms
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
				



   output_file=fopen(OUTPUT,"w");

    if( output_file == NULL )
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }
   fprintf(output_file, "%d\n",N);
   fprintf(output_file, "comment\n");

				for(int ii=0;ii<n*n*n;ii++) printf("atom%d (%f, %f, %f)\n",ii,r[ii][0],r[ii][1],r[ii][2]);
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

      
      for(int ii=0;ii<N;ii++){
         addArray(averageP,averageP,p[ii],3);
         multiplyArray(averageP,averageP,((double)1/N),3);
      }

   }
   saveMomentum("output2",p,N);
   printf("averageP= %f %f %f\n", averageP[0],averageP[1],averageP[2]);

   for(int ss=0;ss<50;ss++){

   for(int jj=0;jj<N;jj++){

   for(int ii=0;ii<N;ii++){

      if(ii!=jj){
      V[ii]+=calculatePotential(r[jj],r[ii],3,EPS,R,L,F);
      calculateForce(force,r[jj],r[ii],3,EPS,R,L,F);
      addArray(sumforce,sumforce,force,3);
      }


   }
   //printf("%f %f %f\n",sumforce[0],sumforce[1],sumforce[2]);
   multiplyArray(sumforce,sumforce,TAU/2,3);
   addArray(p[jj],p[jj],sumforce,3);
   //printf("%f %f %f\n",p[jj][0],p[jj][1],p[jj][2]);
   memset(sumforce, 0, sizeof(sumforce));
   memset(force, 0, sizeof(force));


}

for(int ii=0;ii<N;ii++){
   multiplyArray(force,p[ii],TAU/2/M,3);
   addArray(r[ii],r[ii],force,3);
   
}
 fprintf(output_file, "%d\n",N);
   fprintf(output_file, "comment\n");
savePositions(output_file,r,N);

}

fclose(output_file);


}
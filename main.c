#include <stdio.h>
#include <stdlib.h>
#include "argon.h"

int n;
char *INPUT,*OUTPUT;
int main(int argc, char *argv[]){
	
	INPUT=argv[1];
	OUTPUT=argv[2];
	printf("INPUT FILE: %s\n", INPUT);
	printf("OUTPUT FILE: %s\n", OUTPUT);

	n=loadData(argv[1]);
   printf("liczba atomow: %d\n",n*n*n);



}
#ifndef _argon_h_
#define _argon_h_
#include <stdio.h>


extern void loadData(char *input_file_name,int intParam[],double doubleParam[]);
extern void saveData(char *output_file_name, double data[][3],int N);
extern void saveMomentum(char *output_file_name,double data[][3], int N);
extern void savePositions(FILE *output_file,double data[][3],int N);
extern double distance(double *array1,double *array2, int n);
extern double pressure(double data[][3], double L, int N);
extern double sum(double data[], int N);

extern double temperature(double data[][3],double mass,double k, int n);
extern double kineticEnergy(double data[][3],double mass, int N);

#endif
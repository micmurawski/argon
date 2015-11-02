#ifndef _argon_h_
#define _argon_h_
#include <stdio.h>

struct Parameters{

	int n;
	long double m;
	long double eps;
	long double R;
	long double f;
	long double L;
	long double A;
	long double T;
	long double TAU;
	int s_0;
	int s_d;
	int s_out;
	int s_xyz;
};

extern void loadData(char *input_file_name,struct Parameters *p);
extern void saveData(char *output_file_name, long double data[][3],int N);
extern long double pressure(long double *px,long double *py,long double *pz, long double L, int N);
extern long double temperature(long double *px,long double *py,long double *pz,long double mass,long double k, int n);
extern long double kineticEnergy(long double *px,long double *py,long double *pz,long double mass, int N);
extern long double sumArray(long double *v,int N);

#endif
#ifndef _argon_h_
#define _argon_h_
#include <stdio.h>

struct Parameters{

	int n;
	double m;
	double eps;
	double R;
	double f;
	double L;
	double A;
	double T;
	double TAU;
	int s_0;
	int s_d;
	int s_out;
	int s_xyz;
};

extern void loadData(char *input_file_name,struct Parameters *p);
extern void saveData(char *output_file_name, double data[][3],int N);
extern double pressure(double *px,double *py,double *pz, double L, int N);
extern double temperature(double *px,double *py,double *pz,double mass,double k, int n);
extern double kineticEnergy(double *px,double *py,double *pz,double mass, int N);
extern double sumArray(double *v,int N);
extern void algorytm2(double *x,double *y,double *z,double *Fx,double *Fy,double *Fz,double *V,double EPS, double R, double f,double L,int N);

#endif
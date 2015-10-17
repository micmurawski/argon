#ifndef _argon_h_
#define _argon_h_


void loadData(char *input_file_name,int intParam[],double doubleParam[]);
void saveData(char *output_file_name, double data[][3],int N);
void addArray(double *result,double *array1, double *array2, int n); 
void multiplyArray(double *result,double *array1,double number, int n);
void saveMomentum(char *output_file_name,double data[][3], int N);
double dotProduct(double *array1,double *array2, int n);
double Norm(double *array1,int n);
void subtractArray(double *result,double *array1,double *array2, int n);
double calculatePotential(double *array1, double *array2, int n, double eps, double R);
void calculateForce(double *result, double *array1, double *array2, int n, double eps, double R);


#endif
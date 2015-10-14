#ifndef _argon_h_
#define _argon_h_


void loadData(char *input_file_name,int *data);
void saveData(char *output_file_name, double data[][3],int N);
void addArray(double *result,double *array1, double *array2, int n);
void multiplyArray(double *result,double *array1,double number, int n);
void saveMomentum(char *output_file_name,double data[][3], int N);
double dotProduct(double *array1,double *array2, int n);
double subtractArray(double *result,double *array1,double *array2, int n);


#endif
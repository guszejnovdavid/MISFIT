//Routines for allocating memory for vectors and matrices

#include <stdlib.h>
#include <iostream>

using namespace std;
/*****************************************************************************/
                         /*  Allocation */
//Allocating vector of long integers, size: m
long *vector_long(long m) {
    long *res = (long *) calloc(m, sizeof (long));
    if (res == NULL){
    	cout<<"Error! Failed memory allocation in" <<"vector_long"<<" function.\n";
	}
    return (res);
}

//Allocating matrix of long integers, size: m x n
long **matrix_long(long m, long n) {
    long **res = (long **) calloc(m, sizeof (long *));
    long i;
    for (i = 0; i < m; i++)
        res[i] = (long *) calloc(n, sizeof (long));
    	if (res[i] == NULL){
    		cout<<"Error! Failed memory allocation in" <<"matrix_long"<<" function.\n";
		}
    return (res);
}

//Allocating vector of doubles, size: m
double *vector_double(long m) {
    double *res = (double *) calloc(m, sizeof (double));
    if (res == NULL){
    	cout<<"Error! Failed memory allocation in" <<"vector_double"<<" function.\n";
	}
    return (res);
}

//Allocating vector of booleans, size: m
bool *vector_bool(long m) {
    bool *res = (bool *) calloc(m, sizeof (bool));
    if (res == NULL){
    	cout<<"Error! Failed memory allocation in" <<"vector_bool"<<" function.\n";
	}
    return (res);
}

//Allocating matrix of doubles, size: m x n
double **matrix_double(long m, long n) {
    double **res = (double **) calloc(m, sizeof (double *));
    long i;
    for (i = 0; i < m; i++)
        res[i] = (double *) calloc(n, sizeof (double));
        if (res[i] == NULL){
    		cout<<"Error! Failed memory allocation in" <<"matrix_double"<<" function.\n";
		}
    return (res);
}

//3D matrix of doubles
double ***matrix3d_double(long m, long n, long o){
       long i,j;
       double ***res = (double ***)calloc(m, sizeof(double**));
       for (i = 0; i< m; i++) {
         res[i] = (double **) calloc(n, sizeof(double*));
         for (j = 0; j < n; j++) {
              res[i][j] = (double *)calloc(o, sizeof(double));
            if (res[i][j] == NULL){
    			cout<<"Error! Failed memory allocation in" <<"matrix3d_double"<<" function.\n";
			}
          }
        }
       return (res);
}
/*****************************************************************************/
                         /*  Freeing memory */
//Freeing matrix of long integers, size: mxn
void free_long_matrix(long** M, long m){
     long i;
     for (i = 0; i < m; i++){
         free(M[i]);
         }
         free(M);   
}
//Freeing matrix of doubles, size: mxn
void free_double_matrix(double** M, long m){
     long i;
     for (i = 0; i < m; i++){
         free(M[i]);
         }
         free(M);   
}

void free_double_matrix3d(double*** M, long m, long n){
       long i,j;
       for (i = 0; i< m; i++) {
         for (j = 0; j < n; j++) {
              free(M[i][j]);
          }
          free(M[i]);
        }
}
/*****************************************************************************/
                         /*  Vector Operations */

/* Adding vector and constant */
void vector_set( double x[], double a, long n_elements )
{
 	 //adds up vector plus constant
 	   long i;	     
 	   for(i=0;i<n_elements;i++){
	   							 x[i]=a;	   
 	   
	   }
	   return;
	   }

/* Adding vector and constant */
void vector_const_add( double res[], double x[], double a, long n_elements )
{
 	 //adds up vector plus constant
 	   long i;	     
 	   for(i=0;i<n_elements;i++){
	   							 res[i]=x[i]+a;	   
 	   
	   }
	   return;
	   }

/* Adding two vectors */
void vector_add( double res[], double x1[], double x2[], long n_elements )
{
 	 //adds up two arrays
 	   long i;	     
 	   for(i=0;i<n_elements;i++){
	   							 res[i]=x1[i]+x2[i];	   
 	   
	   }
	   return;
	   }
	   
/*Multiplying vector with constant*/
void vector_const_multiply( double res[], double x[], double factor, long n_elements )
{
 	 //multiplies array with constant
 	   long i;	     
 	   for(i=0;i<n_elements;i++){
	   							 res[i]=x[i]*factor;	   
 	   
	   }
	   return;
	   }
	   
/* Vector dot product */
void vector_dotproduct( double res, double x1[], double x2[], long n_elements )
{
 	 //dot product of two vectors
 	   long i;
       res=0;	     
 	   for(i=0;i<n_elements;i++){
	   							 res+=x1[i]*x2[i];	   
 	   
	   }
	   return;
	   }
/*Multiplying vector with each other for each point*/
void vector_vector_multiply( double res[], double x1[], double x2[], long n_elements )
{
 	 //multiplies array with constant
 	   long i;	     
 	   for(i=0;i<n_elements;i++){
	   							 res[i]=x1[i]*x2[i];	   
 	   
	   }
	   return;
	   }
/*****************************************************************************/
                         /*  Matrix Operations */
                         
/* Adding matrix and vector */
void matrix_const_add( double** res, double** M, double a, long n, long m )
{
 	   long i,j;	     
 	   for(i=0;i<n;i++){
                      for(j=0;j<m;j++){
	   							         res[i][j]=M[i][j]+a;
                                         }
                    }	   
 	   
	   return;
	   }

/* Adding two matrices */
void matrix_add( double** res, double** M1, double** M2, long n, long m )
{
 	   long i,j;	     
 	   for(i=0;i<n;i++){
                      for(j=0;j<m;j++){
	   							         res[i][j]=M1[i][j]+M2[i][j];
                                         }
                    }	   
 	   
	   return;
	   }
	   
/*Multiplying vector with constant*/
void matrix_const_multiply( double** res, double** M, double factor, long n, long m)
{
 	 //multiplies matrix with constant
 	   long i,j;	     
 	   for(i=0;i<n;i++){
                      for(j=0;j<m;j++){
                                       res[i][j]=M[i][j]*factor;
                                       }
 	   
	   }
	   return;
	   }
	   

/*Multiplying matrices with each other for each point*/
void matrix_matrix_multiply(  double** res, double** M1, double** M2, long n, long m )
{
 	 //multiplies array with constant
 	   long i,j;	     
 	   for(i=0;i<n;i++){
                               for(j=0;j<m;j++){
	   							                                       res[i][j]=M1[i][j]*M2[i][j];	   
                                                         }
 	   
	   }
	   return;
	   }

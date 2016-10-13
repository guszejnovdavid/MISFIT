#ifndef vectorfunc_INCLUDED
#define vectorfunc_INCLUDED

/*****************************************************************************/
                         /*  Allocation */
//Allocating vector of long integers, size: m
long *vector_long(long m);

//Allocating matrix of long integers, size: m x n
long **matrix_long(long m, long n);

//Allocating vector of doubles, size: m
double *vector_double(long m);

//Allocating vector of booleans, size: m
bool *vector_bool(long m);

//Allocating matrix of doubles, size: m x n
double **matrix_double(long m, long n);

//Allocating 3D matrix of doubles, size: m x n x o
double ***matrix3d_double(long m, long n, long o);

/*****************************************************************************/
                         /*  Freeing memory */
//Freeing matrix of long integers, size: mxn
void free_long_matrix(long** M, long m);

//Freeing matrix of doubles, size: mxn
void free_double_matrix(double** M, long m);

//Freeing 3D matrix of doubles, size: m x n x o
void free_double_matrix3d(double*** M, long m, long n);

/*****************************************************************************/
                         /*  Vector Operations */

//Setting vector to constant
void vector_set( double x[], double a, long n_elements );

/* Adding vector and constant */
void vector_const_add( double res[], double x[], double a, long n_elements );

/* Adding two vectors */
void vector_add( double res[], double x1[], double x2[], long n_elements );

/*Multiplying vector with constant*/
void vector_const_multiply( double res[], double x[], double factor, long n_elements );
   
/* Vector dot product */
void vector_dotproduct( double res, double x1[], double x2[], long n_elements );

/*Multiplying vector with each other for each point*/
void vector_vector_multiply( double res[], double x1[], double x2[], long n_elements );

/*****************************************************************************/
                         /*  Matrix Operations */
                         
/* Adding matrix and constant */
void matrix_const_add( double** res, double** M, double a, long n, long m );

/* Adding two matrices */
void matrix_add( double** res, double** M1, double** M2, long n, long m );
	   
/*Multiplying vector with constant*/
void matrix_const_multiply( double** res, double** M, double factor, long n, long m);

/*Multiplying matrices with each other for each point*/
void matrix_matrix_multiply(  double** res, double** M1, double** M2, long n, long m );


#endif

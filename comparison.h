#ifndef comparison_INCLUDED
#define comparison_INCLUDED
// Comparing variables

// Minimum of 2 doubles
double Min_2_double(double a, double b);

 
// Maximum of 2 doubles
double Max_2_double(double a, double b);

// Minimum of 2 longs
double Min_2_long(long a, long b);
 
// Maximum of 2 longs
double Max_2_long(long a, long b);

// Is a within [xmin,xmax]? If yes return a else give the limits.
double InLimits_double(double a, double xmin, double xmax);

// Minimum of array of doubles
double Min_array_double(double* x, long N);

// Maximum of array of doubles
double Max_array_double(double* x, long N);
#endif

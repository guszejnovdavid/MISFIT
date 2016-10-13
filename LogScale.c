//Create vector with logarithmic spacing

#include <math.h>
#include "LogScale.h"

void LogScale(double* vector, double x0, double x1, long n)
{
long i; //loop variable
double step=log(x1/x0)/(n-1);
double logx0=log(x0);

for (i = 0; i < n; i++){
    vector[i]=exp(logx0+i*step);
    }
}

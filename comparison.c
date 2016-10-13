// Collection of routines for comparing variables

#include "comparison.h"

// Minimum of 2 doubles
double Min_2_double(double a, double b)
{
    if (b<a){
       return b;
       }
       else
       {
       return a;
       }
}
 
// Maximum of 2 doubles
double Max_2_double(double a, double b)
{
    if (b<a){
       return a;
       }
       else
       {
       return b;
       }
}

// Minimum of 2 longs
double Min_2_long(long a, long b)
{
    if (b<a){
       return b;
       }
       else
       {
       return a;
       }
}
 
// Maximum of 2 longs
double Max_2_long(long a, long b)
{
    if (b<a){
       return a;
       }
       else
       {
       return b;
       }
}

// Is a within [xmin,xmax]? If yes return a else give the limits.
double InLimits_double(double a, double xmin, double xmax)
{
    if (a>xmax){
       return xmax;
       }
       else
       {
         if(a<xmin){
                return xmin;
                }
         else
         {
           return a;
          }
       }
}

// Minimum of array of doubles
double Min_array_double(double* x, long N)
{
long i;
double res;

res=x[0]; //init
 for (i = 0; i < N; i++){
     if (x[i]>res){
        res=x[i];
        }
     }
 return res;
}

// Maximum of array of doubles
double Max_array_double(double* x, long N)
{
long i;
double res;

res=x[0]; //init
 for (i = 0; i < N; i++){
     if (x[i]<res){
        res=x[i];
        }
     }
 return res;
}

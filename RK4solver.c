#include <stdlib.h>
#include "vectorfunctions.h"
#include "Deriv_func.h"
#include "RK4solver.h"

void rk4_step ( double y[], double t, double dt, long dimension )


/******************************************************************************/
/*
  Purpose:

    Carries out the 4th order Runge-Kutta integration between t and tout
*/
{

  long j; //loop integer
  double *temp, *k1, *k2, *k3, *k4; //intermediate results
  //Allocating memory for arrays
  temp = ( double * ) malloc ( dimension * sizeof ( double ) );
  k1 = ( double * ) malloc (dimension * sizeof ( double ) );
  k2 = ( double * ) malloc ( dimension * sizeof ( double ) );
  k3 = ( double * ) malloc ( dimension * sizeof ( double ) );
  k4 = ( double * ) malloc ( dimension * sizeof ( double ) );
  
 
 Deriv_func ( t , y, k1, dt, dimension ); //calculating k1
 
 vector_const_multiply( k1, k1, dt, dimension );
 vector_const_multiply( temp, k1, 0.5, dimension ); // getting 0.5*k1
 vector_add( temp, y, temp,dimension ); // getting y+0.5*k1
 
 Deriv_func ( (t+0.5*dt) , temp, k2, dt, dimension ); //calculating k2
 
 vector_const_multiply( k2, k2, dt,dimension ); 
 vector_const_multiply( temp, k2, 0.5,dimension ); // getting 0.5*k2
 vector_add( temp, y, temp,dimension ); // getting y+0.5*k2
 Deriv_func ( (t+0.5*dt) , temp, k3, dt, dimension ); //calculating k3
 
 vector_const_multiply( k3, k3, dt,dimension ); 
 vector_add( temp, y, k3,dimension ); // getting y+k3
 Deriv_func ( (t+dt) , temp, k4, dt, dimension ); //calculating k4
 vector_const_multiply( k4, k4, dt,dimension );
 //Getting next y values
 for(j=0;j<dimension;j++){
 					 y[j]+=(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6.0;
					 }		  

//Cleanup
free(temp);
free(k1);
free(k2);
free(k3);
free(k4);

return;           

      }

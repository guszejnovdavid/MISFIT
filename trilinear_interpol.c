//Trilinear interpolation
/*
Copied from the source code by A. Odrzywolek
*/

#include <math.h>
#include "trilinear_interpol.h"



inline double trilinear_interp_unit_cube
(double X, double Y, double Z, double A000, double A100, double A101, double A001, double A010, double A110, double A111, double A011)
{

  return A000*(1.0-X)*(1.0-Y)*(1.0-Z) + A100*X*(1.0-Y)*(1.0-Z) + A010*(1.0-X)*Y*(1.0-Z) + A001*(1.0-X)*(1.0-Y)*Z 
       + A101* X*(1.0 - Y)*Z + A011* (1.0 - X)*Y*Z + A110* X* Y* (1.0 - Z) + A111* X*Y*Z;
}

double trilinear_interp(double X, double Y, double Z, const double *tbl, 
      const int N_X, const int N_Y, const int N_Z,
      const double X_MIN, const double X_MAX, 
      const double Y_MIN, const double Y_MAX,
      const double Z_MIN, const double Z_MAX,
      const int OUT_OF_DOMAIN, const double OUT_OF_DOMAIN_VALUE)
{

  double f111,f211,f212,f112, f121,f221,f222,f122;
  double X_A,X_B,Y_A,Y_B,Z_A,Z_B;

/* grid spacing */
  const double delta_X = (X_MAX-X_MIN)/(N_X-1);
  const double delta_Y = (Y_MAX-Y_MIN)/(N_Y-1);
  const double delta_Z = (Z_MAX-Z_MIN)/(N_Z-1);
 
  const double delta_X_inv = (N_X-1)/(X_MAX-X_MIN);
  const double delta_Y_inv = (N_Y-1)/(Y_MAX-Y_MIN);
  const double delta_Z_inv = (N_Z-1)/(Z_MAX-Z_MIN);
  int i,j,k;
  double X_unit, Y_unit, Z_unit;


  /* select upper right grid point BUG!!!!!! randomly (but very rarely) 
 sometimes select WRONG point!!!*/
  i = (int) ceil( (X-X_MIN)*delta_X_inv);
  j = (int) ceil( (Y-Y_MIN)*delta_Y_inv);
  k = (int) ceil( (Z-Z_MIN)*delta_Z_inv);

/* Return OUT_OF_DOMAIN_VALUE if X or Y or Z are under/overflown */
/* edge of the domain and out-of-range values */
   if( OUT_OF_DOMAIN==0 ) 
   {
     if( (X<X_MIN) || (Y<Y_MIN) || (Z<Z_MIN) || (X>X_MAX) || (Y>Y_MAX) || (Z>Z_MAX) ) return OUT_OF_DOMAIN_VALUE;
   };
   
/* cell corner values */ 
// f222=TBL[i][j][k] declared as double TBL[N_X+1][N_Y+1][N_Z+1]; tbl=&TBL[0][0][0]

  f111 = *(tbl+(N_Z)*(N_Y)*(i-1)+(N_Z)*(j-1)+k-1)   ; //TBL[i-1][j-1][k-1];
  f211 = *(tbl+(N_Z)*(N_Y)*i+(N_Z)*(j-1)+k-1)   ;// TBL[i][j-1][k-1];
  f212 = *(tbl+(N_Z)*(N_Y)*i+(N_Z)*(j-1)+k)   ;//TBL[i][j-1][k];
  f112 = *(tbl+(N_Z)*(N_Y)*(i-1)+(N_Z)*(j-1)+k)   ;//TBL[i-1][j-1][k];
  f121 = *(tbl+(N_Z)*(N_Y)*(i-1)+(N_Z)*j+k-1);    //TBL[i-1][j][k-1];
  f221 = *(tbl+(N_Z)*(N_Y)*i+(N_Z)*j+k-1);  //TBL[i][j][k-1];
  f222 = *(tbl+(N_Z)*(N_Y)*i+(N_Z)*j+k)   ; //TBL[i][j][k];
  f122 = *(tbl+(N_Z)*(N_Y)*(i-1)+(N_Z)*j+k); //TBL[i-1][j][k];

  X_A = X_MIN + (i-1)*delta_X;    
  X_B = X_MIN + i*delta_X;    
  
  Y_A = Y_MIN + (j-1)*delta_Y;
  Y_B = Y_MIN + j*delta_Y;

  Z_A = Z_MIN + (k-1)*delta_Z;
  Z_B = Z_MIN + k*delta_Z;


/* Rescalling to unit cube {0,1}^3 */
  X_unit = (X-X_A)/(X_B-X_A);
  Y_unit = (Y-Y_A)/(Y_B-Y_A);
  Z_unit = (Z-Z_A)/(Z_B-Z_A);

  return trilinear_interp_unit_cube(X_unit,Y_unit,Z_unit, f111,f211,f212,f112,f121,f221,f222,f122);

}

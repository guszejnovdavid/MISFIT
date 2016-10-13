// Newton method for finding Mach number. Algorith modified from code by John H. Mathews 1995


/*
---------------------------------------------------------------------------
Algo2-5.c   C program for implementing Algorithm 2.5
Algorithm translated to  C  by: Dr. Norman Fahrer
IBM and Macintosh verification by: Daniel Mathews

NUMERICAL METHODS: C Programs, (c) John H. Mathews 1995
To accompany the text:
NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
Prentice Hall, Inc.; USA, Canada, Mexico ISBN 0-13-624990-6
Prentice Hall, International Editions:   ISBN 0-13-625047-5
This free software is compliments of the author.
E-mail address:       "mathews@fullerton.edu"

Algorithm 2.5 (Newton-Raphson Iteration).
Section   2.4, Newton-Raphson and Secant Methods, Page 84
---------------------------------------------------------------------------
*/
/*
---------------------------------------------------------------------------

Algorithm 2.5 (Newton-Raphson Iteration). To find a root
f(x) = 0 given one initial approximation  p_0  and using the iteration

                                                  f(p_(k-1))  
      p_k  =  p_(k-1)  -  -------------     for k = 1, 2, ...  
                                                 f'(p_(k-1))  

---------------------------------------------------------------------------
*/

/* User has to supply a function named :  ffunction
   and its first derivative :            dffunction

   An example is included in this program */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "FindParameter_Newton.h"
//Include global variables
#include "GlobalVars.h"

/******************************************************************************/
							/*  Mass-Mach number relation   */

    double ffunction_Mach(double x, double m)
    {
        return (Ms*x*x*(1+x*x)/(1+1.0/x)-m); //Mach number-mass relation
    }

/*  differential of Mass-Mach number relation */

    double dffunction_Mach(double x)
    {
        //return ( 3.0 * ( (3.0*x*x+5.0*x*x*x*x)/(1.0+x)-(x*x*x)*(1+x*x)/(1.0+x)/(1.0+x) ) ); //Mach number-mass relation
        return (  (Ms*x*x*(3.0 + x*(2.0 + x*(5.0 + 4.0*x))))/((1 + x)*(1+x)) );
    }


/* -------------------------------------------------------- */

/*  Find Mach number  */

    double FindMach(double m)

{
    double Delta = 1E-6;       /* Tolerance  */
    double Epsilon = 1E-6;     /* Tolerance  */
    double Small = 1E-6;       /* Tolerance  */

    int Max = 99;  /* Maximum number of iterations  */
    int Cond = 0;  /* Condition fo loop termination */
    int K;         /* Counter for loop              */

    double P0=1.0;    /* INPUT : Must be close to the root */
    double P1;    /* New iterate    */
    double Y0;    /* Function value */
    double Y1;    /* Function value */
    double Df;    /* Derivative     */
    double Dp;
    double RelErr;


    Y0 = ffunction_Mach(P0,m);

    for ( K = 1; K <= Max ; K++) {

        if(Cond) break;

        Df = dffunction_Mach(P0);      /* Compute the derivative  */

         if( Df == 0) {       /* Check division by zero */
             Cond = 1;
             Dp   = 0;
         }

         else Dp = Y0/Df;

        P1 = P0 - Dp;       /* New iterate */
        Y1 = ffunction_Mach(P1,m);   /* New function value */

        RelErr = 2 * fabs(Dp) / ( fabs(P1) + Small );  /* Relative error */

        if( (RelErr < Delta)  && (fabs(Y1) < Epsilon) ) { /* Check for   */

            if( Cond != 1) Cond = 2;                      /* convergence */

        }

        P0 = P1;
        Y0 = Y1;
    }

if(Cond == 0) printf("The maximum number of iterations was exceeded !\n");

if(Cond == 1) printf("Division by zero was encountered !\n");

return P0;

//if(Cond == 2) printf("The root was found with the desired tolerance !\n");

}   /* End of main program */


/******************************************************************************/
							/*  Mass-Radius relation   */

    double ffunction_R(double x, double m)

    {
        return (Ms/2.0*(x/Rs)*(1.0+x/Rs)-m); //R-mass relation
    }

/*  differential of Mass-Mach number relation */

    double dffunction_R(double x)

    {
        return ( Ms/(2.0*Rs)*( (1.0+x/Rs) + (x/Rs)  ) ); //R-mass relation derivative
    }


/* -------------------------------------------------------- */

/*  Find Mach number  */

    double FindR(double m)

{
    double Delta = 1E-6;       /* Tolerance  */
    double Epsilon = 1E-6;     /* Tolerance  */
    double Small = 1E-6;       /* Tolerance  */

    int Max = 99;  /* Maximum number of iterations  */
    int Cond = 0;  /* Condition fo loop termination */
    int K;         /* Counter for loop              */

    double P0=1.0;    /* INPUT : Must be close to the root */
    double P1;    /* New iterate    */
    double Y0;    /* Function value */
    double Y1;    /* Function value */
    double Df;    /* Derivative     */
    double Dp;
    double RelErr;


    Y0 = ffunction_R(P0,m);

    for ( K = 1; K <= Max ; K++) {

        if(Cond) break;

        Df = dffunction_R(P0);      /* Compute the derivative  */

         if( Df == 0) {       /* Check division by zero */
             Cond = 1;
             Dp   = 0;
         }

         else Dp = Y0/Df;

        P1 = P0 - Dp;       /* New iterate */
        Y1 = ffunction_R(P1,m);   /* New function value */

        RelErr = 2 * fabs(Dp) / ( fabs(P1) + Small );  /* Relative error */

        if( (RelErr < Delta)  && (fabs(Y1) < Epsilon) ) { /* Check for   */

            if( Cond != 1) Cond = 2;                      /* convergence */

        }

        P0 = P1;
        Y0 = Y1;
    }

if(Cond == 0) printf("The maximum number of iterations was exceeded !\n");

if(Cond == 1) printf("Division by zero was encountered !\n");

return P0;

//if(Cond == 2) printf("The root was found with the desired tolerance !\n");

}   /* End of main program */




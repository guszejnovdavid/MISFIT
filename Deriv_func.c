#include <math.h>
#include "gammafunc.h"
#include "Deriv_func.h"
#include "GlobalVars.h"

void Deriv_func ( double t, double y[], double yp[], double dt, long dimension )

/******************************************************************************/
/*
  Purpose:
  		  
    Evaluates the derivative, using equation for scalar size
*/
{
   /*Calculating derivative*/ 
   // In this case y[0]= r relative radius, y[1]=sigma surface density, y[2]=edge Mach number at the beginning, y[3] is the current edge Mach number squared
      // d r/d tau=1/(1-Q/2) -1/sqrt(r) Mach^3/(1+Mach^2)^(3/2)*sqrt(M/M_0), see http://arxiv.org/abs/1507.06678 (last correction factor accounts for mass loss)
   //Safegurad against negative values which cause problems
   if(y[0]>0){
   	yp[0]=-1.0/(1.0-QVirial/2.0)*pow(y[0],-0.5)*pow(y[2]/(1.0+y[2]),1.5)*sqrt(y[3]);
      // Sigma, evolved separetely
    yp[1]=0.0;
     // Edge Mach number, evolved separetely
   	yp[2]=0.0;
     // Mass ratio, evolved separetely
   	yp[3]=0.0;
   }
   else{
   	yp[0]=0;
   	yp[1]=0;
   	yp[2]=0;
   	yp[3]=0;
   }
   
   
  return;
}

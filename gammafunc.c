/**************************************************************************************/
								/*  MISFIT Equation of State */
								
// This contains the equation of state for the medium in the form gamma=gamma(...)

#include <math.h>
#include "gammafunc.h"
#include "GlobalVars.h"

#define Rhocrit1 5E3 //critical density for sub-isothermal->isothermal transition for the non scaled EOS in msolar/pc^3
#define Rhocrit2 5E8 //critical density for isothermal->adiabatic transition for the non scaled EOS in msolar/pc^3

double gammafunc (double rhoavg, double T )

//Calculates polytropic index for a given volume density
{
   double res;

#if (IsothermalEOS==false)   
	if ( rhoavg<Rhocrit1 ) 
	{
	    res=0.7;// thin, sub-isothermal regime
	}
	else 
	{
		if ( rhoavg<Rhocrit2 ) // lower limit
		{
		    res=1.0;// isothermal regime
		}
		else{
			res=1.4; // opaque, adiabatic regime
		}
	}
#else
	//Force isothermal collapse
	res=1.0;
#endif

  return res;
}

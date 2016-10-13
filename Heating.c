/**************************************************************************************/
								/*  MISFIT Protostellar Heating */
								
// Manages the calculation (rank=0 process is the Master), collects the results and distributes tasks.

#include <math.h>
#include "Heating.h"
#include "GlobalVars.h"

#define TempPrefactor 0.595403 // T_heat=TempPrefactor*M^(3/8)*R^(-7/8)


double HeatingTemp ( double Mcloud, double Rcloud, double heatingeff )

/******************************************************************************/
/*
  Purpose:
  		  
    Calculates the temperature to which the protostar seed at the center of the cloud would heat up the gas (assuming homogenity)
    Formula: T^4=Psi*Sqrt(G)*M^(3/2)/[4*pi*sigma*R^(7/2)]
*/
{  
  #if (HeatingAllowed==true)
  return heatingeff*TempPrefactor*pow(Mcloud,0.375)*pow(Rcloud,-0.875);
  #else
  return 0;
  #endif
}

void RhoTempFieldWithHeating(double HeatTemp, double* TField, double* RhoField, double* Rho_T_Field, long N)
// Calculates the rho*T field where heating is taken into account
{
	long i;
	double factor=1.0;
	double HeatTemp4=pow(HeatTemp,4.0);
	#if (R05Heating==true)
	double x,y,z;
	const long Nsq=Ngrid*Ngrid;
	const double Nhalf=(double)(Ngrid-1)/2.0;
	#endif
	
	for(i=0;i<N;i++){
		#if (R05Heating==true)
		x= (double)(i/Nsq)-Nhalf; //x coordinate index
     	y= (double)((i % Nsq)/Ngrid)-Nhalf; //y coordinate index
     	z= (double)(i % Ngrid)-Nhalf; //z coordinate index
     	factor=Nhalf*Nhalf/(x*x+y*y+z*z); // (sqrt(r0/r))^4
     	#endif
		Rho_T_Field[i]=pow(factor*HeatTemp4+pow(TField[i],4.0),0.25)*RhoField[i];
	}
}

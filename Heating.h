#ifndef Heating_INCLUDED
#define Heating_INCLUDED

double HeatingTemp ( double Mcloud, double Rcloud, double heatingeff  );
void RhoTempFieldWithHeating(double HeatTemp, double* TField, double* RhoField, double* Rho_T_Field, long N);

#endif

// Newton method for finding Mach number. Algorith modified from code by John H. Mathews 1995

#ifndef FindParameter_INCLUDED
#define FindParameter_INCLUDED
// Mach-mass relation
    double ffunction_Mach(double x, double m);
//Derivative
    double dffunction_Mach(double x);
// Mach finder
    double FindMach(double m);
// R-mass relation
    double ffunction_R(double x, double m);
//Derivative
    double dffunction_R(double x);
// R finder
    double FindR(double m);
#endif

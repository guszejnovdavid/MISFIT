#ifndef EvolveField_INCLUDED
#define EvolveField_INCLUDED

void EvolveField(struct inputstruct_s* inputstruct);

double SoundSpeedSquare(double TMean, double gamma);

 //Critical density using Jeans criterion     
inline double RhoCrit(double TMean, double R, double MachedgeSq, double gamma);

inline bool IsSelfGravitating(double SphereMass, double TMean, double R, double MachedgeSq, double gamma, double StrictnessFactor);

double FindCMFTime(double tcurrent);

long FindCMFTimeIndex(double tcurrent);

double SetOptimaldt(double dt,double renormstepMax,double dDeltaKMax);

#endif

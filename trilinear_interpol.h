//Trilinear interpolation

#ifndef trilinear_INCLUDED
#define trilinear_INCLUDED

inline double trilinear_interp_unit_cube
(double X, double Y, double Z, double A000, double A100, double A101, double A001, double A010, double A110, double A111, double A011);

double trilinear_interp(double X, double Y, double Z, const double *tbl, 
      const int N_X, const int N_Y, const int N_Z,
      const double X_MIN, const double X_MAX, 
      const double Y_MIN, const double Y_MAX,
      const double Z_MIN, const double Z_MAX,
      const int OUT_OF_DOMAIN, const double OUT_OF_DOMAIN_VALUE);

#endif

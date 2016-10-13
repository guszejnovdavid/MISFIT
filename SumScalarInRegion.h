/*
Calculates the sum of a scalar field in different 3D regions. The background must be cubic;


*/

#ifndef SumScalar_INCLUDED
#define SumScalar_INCLUDED

// Sum in a sphere
double SumScalar_Sphere(
       double* Field, long N, //field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, // limits of grid
       double R, double x0, double y0, double z0); //parameters of sphere

//Average in a speher
double AverageScalar_Sphere(
       double* Field, long N, //field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, // limits of grid
       double R, double x0, double y0, double z0); //parameters of sphere

// Sum in a cube
double SumScalar_Slab(
       double* Field, long N,//field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax,// limits of grid
       double a, double x0, double y0, double z0); //parameters of cube
       
// Average in a cube
double AverageScalar_Slab(
       double* Field, long N,//field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax,// limits of grid
       double a, double x0, double y0, double z0); //parameters of cube

#endif

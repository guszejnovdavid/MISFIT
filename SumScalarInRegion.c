#include "SumScalarInRegion.h"
#include "comparison.h"
#include <math.h>
/*
Calculates the sum of a scalar field in different 3D regions. The background must be cubic;


*/

// Sum in a sphere
double SumScalar_Sphere(
       double* Field, long N, //field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, // limits of grid
       double R, double x0, double y0, double z0) //parameters of sphere
{
double res,temp;
long i,j,k;
long Imin,Imax,Jmin,Jmax,Kmin,Kmax; //indices
double xsq_norm, ysq_norm/*, zsq_norm*/ ; //used in distance calculations
double rsq; //radial distance square of point to origin normalized to d
double contribfactor; //fraction of volume element within the sphere, <1.0 for boundary points
double reducedr; // (R-(r-d/2))/d, used for getting contribfactor

const long Nsq=N*N; //N^2 used for indexing
const double d=(Xmax-Xmin)/(double)(N-1); //spacing
//const double Rsq=R*R/(d*d); //radius squared renormalized to grid
const double dV=d*d*d; //volume element, assuming cubic lattice
const double Rprime=R+d/2.0; //actual radius within which points contribute
const double Rprimesq=Rprime*Rprime/(d*d); //modified radius squared renormalized to grid
const double Rfull=R-d/2.0; //radius within which points contribute at 100%
const double Rfullsq=Rfull*Rfull/(d*d); //modified radius squared renormalized to grid
// Origin coordinates
const double x0_norm=(x0-Xmin)/d;
const double y0_norm=(y0-Ymin)/d;
const double z0_norm=(z0-Zmin)/d;

res=0; //final result

//Get smaller cube where the points can be
  //X index
Imin=(long)Max_2_double(floor((InLimits_double(x0-Rprime, Xmin, Xmax)-Xmin)/d), 0.0);
Imax=(long)Min_2_double(ceil((InLimits_double(x0+Rprime, Xmin, Xmax)-Xmin)/d),(double)(N-1));

 //Go through all points in subcube
for (i = Imin; i <=Imax; i++){
    xsq_norm=(i-x0_norm)*(i-x0_norm);
    temp=sqrt(Rprimesq-xsq_norm);
	if(temp>0){ //check to make sure we don't get negative values
	//Y index
	Jmin=(long)Max_2_double(floor(y0_norm-temp), 0.0);
	Jmax=(long)Min_2_double(ceil(y0_norm+temp), (double)(N-1));
    for (j = Jmin; j <=Jmax; j++){
        ysq_norm=(j-y0_norm)*(j-y0_norm);
        temp=sqrt(Rprimesq-xsq_norm-ysq_norm);
		if(temp>0){ //check to make sure we don't get negative values
		//Z index
		Kmin=(long)Max_2_double(floor(z0_norm-temp), 0.0);
		Kmax=(long)Min_2_double(ceil(z0_norm+temp), (double)(N-1));
        for (k = Kmin; k <=Kmax; k++){
        	//zsq_norm=(k-z0_norm)*(k-z0_norm);
        	rsq=xsq_norm+ysq_norm+(k-z0_norm)*(k-z0_norm);
        	if(Rfullsq>rsq){ //close enough to contribute 100%
        		res+=Field[k+j*N+i*Nsq];
			}
			else{
				if(Rprimesq>rsq){ //close enough to contribute but not 100%
					reducedr=R/d-sqrt(rsq)+0.5; // (R-(r-d/2))/d used for getting contribfactor, between 0 and 1
					contribfactor=pow(reducedr,3.0)*(10.0-15.0*reducedr+6.0*reducedr*reducedr); //10x^3+6x^5-15x^4
        			res+=contribfactor*Field[k+j*N+i*Nsq];
				}	
			}
         }
     }//end if
     }
 }//end if
}

              
 return (res*dV); //return sum of scalar field             
              
}

//Average in a sphere
double AverageScalar_Sphere(
       double* Field, long N, //field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, // limits of grid
       double R, double x0, double y0, double z0) //parameters of sphere
{
double res;
long i,j,k;
long Counted=0;
long Imin,Imax,Jmin,Jmax,Kmin,Kmax; //indices
double xsq_norm, ysq_norm, zsq_norm; //used in distance calculations
double temp; //temp number

const long Nsq=N*N; //N^2 used for indexing
const double d=(Xmax-Xmin)/(double)(N-1); //spacing
const double Rsq=R*R/(d*d); //radius squared renormalized to grid
// Origin coordinates
const double x0_norm=(x0-Xmin)/d;
const double y0_norm=(y0-Ymin)/d;
const double z0_norm=(z0-Zmin)/d;

res=0; //final result


//Get smaller cube where the points can be
  //X index
Imin=lround((InLimits_double(x0-R, Xmin, Xmax)-Xmin)/d);
Imax=lround((InLimits_double(x0+R, Xmin, Xmax)-Xmin)/d);
              
// Go through all points in subcube
for (i = Imin; i < (Imax+1); i++){
    xsq_norm=(i-x0_norm)*(i-x0_norm);
    temp=sqrt(Rsq-xsq_norm);
	if(temp>0.0){ //check to make sure we don't get negative values
	//Y index
	Jmin=(long)ceil(y0_norm-temp);
	Jmax=(long)floor(y0_norm+temp);
    for (j = Jmin; j < (Jmax+1); j++){
        ysq_norm=(j-y0_norm)*(j-y0_norm);
        temp=sqrt(Rsq-xsq_norm-ysq_norm);
		if(temp>0.0){ //check to make sure we don't get negative values
		//Z index
		Kmin=(long)ceil(z0_norm-temp);
		Kmax=(long)floor(z0_norm+temp);
        for (k = Kmin; k < (Kmax+1); k++){
        	zsq_norm=(k-z0_norm)*(k-z0_norm);
        	if(Rsq>xsq_norm+ysq_norm+zsq_norm){
              				 ++Counted;
                             res+=Field[k+j*N+i*Nsq];
            }
         }
     }//end if
     }
 }//end if
}

if(Counted){        
 return (res/Counted); //return average of scalar field 
}
else{
	return 0;
}
              
}


// Sum in a cube
double SumScalar_Slab(
       double* Field, long N,//field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax,// limits of grid
       double a, double x0, double y0, double z0) //parameters of cube
{
double res;
long i,j,k;
long Imin,Imax,Jmin,Jmax,Kmin,Kmax; //indices

const long Nsq=N*N; //N^2 used for indexing
const double d=(Xmax-Xmin)/(double)(N-1); //spacing
const double dV=d*d*d; //volume element, assuming cubic lattice

res=0; //final result


//Get smaller cube where the points can be
  //X index
Imin=lround((InLimits_double(x0-a/2.0, Xmin, Xmax)-Xmin)/d);
Imax=lround((InLimits_double(x0+a/2.0, Xmin, Xmax)-Xmin)/d);
  //Y index
Jmin=lround((InLimits_double(y0-a/2.0, Ymin, Ymax)-Ymin)/d);
Jmax=lround((InLimits_double(y0+a/2.0, Ymin, Ymax)-Ymin)/d);
  //Z index
Kmin=lround((InLimits_double(z0-a/2.0, Zmin, Zmax)-Zmin)/d);
Kmax=lround((InLimits_double(z0+a/2.0, Zmin, Zmax)-Zmin)/d);
              
// Go through all points in subcube
for (i = Imin; i < (Imax+1); i++){
    for (j = Jmin; j < (Jmax+1); j++){
        for (k = Kmin; k < (Kmax+1); k++){
             res+=Field[k+j*N+i*Nsq]*dV;
         }
     }
}
              
 return res; //return sum of scalar field             
              
}

// Average in a cube
double AverageScalar_Slab(
       double* Field, long N,//field array and length of one dimension
       double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax,// limits of grid
       double a, double x0, double y0, double z0) //parameters of cube
{
double res;
long i,j,k;
long Counted=0;
long Imin,Imax,Jmin,Jmax,Kmin,Kmax; //indices

const long Nsq=N*N; //N^2 used for indexing
const double d=(Xmax-Xmin)/(double)(N-1); //spacing

res=0; //final result


//Get smaller cube where the points can be
  //X index
Imin=lround((InLimits_double(x0-a/2.0, Xmin, Xmax)-Xmin)/d);
Imax=lround((InLimits_double(x0+a/2.0, Xmin, Xmax)-Xmin)/d);
  //Y index
Jmin=lround((InLimits_double(y0-a/2.0, Ymin, Ymax)-Ymin)/d);
Jmax=lround((InLimits_double(y0+a/2.0, Ymin, Ymax)-Ymin)/d);
  //Z index
Kmin=lround((InLimits_double(z0-a/2.0, Zmin, Zmax)-Zmin)/d);
Kmax=lround((InLimits_double(z0+a/2.0, Zmin, Zmax)-Zmin)/d);
              
// Go through all points in subcube
for (i = Imin; i < (Imax+1); i++){
    for (j = Jmin; j < (Jmax+1); j++){
        for (k = Kmin; k < (Kmax+1); k++){
        	 Counted++;
             res+=Field[k+j*N+i*Nsq];
         }
     }
}
              
if(Counted){        
 return (res/(double)Counted); //return average of scalar field 
}
else{
	return 0;
}           
              
}

/**************************************************************************************/
			/*  Main Evolution Routine for MISFIT */
			
// This evolves the density and temperature fields of a collapsing gas cloud while looking for self-gravitating structures.
// If such a fragment is found, it is removed and sent to the Master which will send it to another slave routine that will evolve it using this routine. 

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>

//Include own stuff
#include "comparison.h"
#include "LogScale.h"
#include "vectorfunctions.h"
#include "EvolveField.h"
#include "RK4solver.h"
#include "SumScalarInRegion.h"
#include "gammafunc.h"
#include "trilinear_interpol.h"
#include "gaussrand.h"
#include "fft.h"
#include "Restart.h"
#include "Heating.h"
//Include global variables
#include "GlobalVars.h"
#include "FileManage.h"
#include "InputDataManage.h"

//Make it run with MPI for cluster
# include "mpi.h"

    
using namespace std;
const double CMFTimes[CMFTIMENUM] = CMFTIMEVALUES; //times at which CMF and IMF are saved.
const double FourPiOver3=4.0/3.0*M_PI; //shorthand for 4*p/3 


//Evolves the T and rho fields of a cloud as it collapses and checks for fragmentation
void EvolveField(struct inputstruct_s* inputstruct)
{
	
//Input variables
double* rho_field; //the density field
double*  T_field; //temperature field
double R0=inputstruct[0].R0; //initial size
double X0, Y0, Z0; //initial position
double Machedge;//initial Mach number
double Vx, Vy, Vz; //initial velocity
double InitScale; //relative scale below which turbulence is initialized, units of R0
double tparent=inputstruct[0].tparent; //time at start
     
//Constants
const long NSphere0=2; //Base number of spheres to use in Mont-Carlo
const double MCFillFactor=2.0; //filling factor, volume of MC spheres/original volume
const long Ngridsq=Ngrid*Ngrid; // shorthand, used for indexing
const double gridspacing=2.0/(Ngrid-1); //Grid spacing in R units
const double Kunit=2*M_PI/Ngrid; //Smallest k-> largest mode
const double KSpaceDensity=pow(Kunit,-3.0); //Density in k space
//const double Qprime=1.0;// virial parameter used in collapse threhold
const double RhoEmpty=1E-10; //density to use for empty space, set as 1E-10 Msolar/pc^3
const double TEmpty=10; //tempearture to use for empty space, set as 10 K

//Cloud contraction parameters
//Parameters
const long MaxStepNum=100001; //maximum number of steps
const double dt0=DEFAULTTIMESTEP; //default time step
double Tmax=10.0; //Maximum time for evolution
const double MinStep=Tmax/(double)(MaxStepNum-1); //smallest step allowed
const double ymin=MINIMUMSIZE; //Relative size of the cloud when calculation stops.
const double MaxAbsErr=R0*MINIMUMSIZE/100.0; //Maximum absolute error allowed during a timestep
const long dimension=4; //dimension of RK4 diff. eq. (#variables+ #parameters)

double CMF_Time; //time at which we are interested in the CMF, using Myr units
long CMF_Time_Index; //index corresponding to CMF_Time

//Helping variables for cloud contraction
double dt=dt0,t=0.0, err; // timestep, current time, error
double tcurrent=tparent; //current time in units where t_dynamical=sqrt(R^3/M), so d t_current= dt*sqrt(R^3/M)
double timeconversionfactor; //shorthand for frefall time, sqrt(R0^3/M0), in Myr
double tcrossing; //crossing time, in Myr
double timescalecorrection=1.0; //used for the ratio: t_freefall(t=0)/t_crossing(t), used in density field mode evolution
double *y0,*y1,*y2; //possible results of RK4 steps
ofstream myfile; //Output file handle

//Output structure
struct inputstruct_s *outputstruct; // to store data for substructure to be used for recursion

#if (DensitySave==true)
	//structure used to send field data to master
	struct densitystruct_s *fieldstruct;
#endif

//MPI communication
int message; //message code 
//MPI_Status status; //MPI message status
int mpierror; //integer for MPI error messages
int rank; // process rank

//Variables
long i,j,k,l; // loop variables
double temp; //temp var
//double *temppointer; //temporary pointer
double R=R0; //Current radius
double Mass, MassInit, MassWithinRadius; //total mass, initial mass, mass within the cloud (counting fragments)
double RhoMean; //mean density
double TMean; //mean temperature
double RhoMean0; //original mean density
double NormFactor; //factor for normalizing density; used in the main loop
double *RelSizeVector; //relative size of cloud
double gamma=1.0; //polytropic index of cloud
double *S; //S, variance for each Fourier scale
double *Mach; //Mach number for turbulent velocity dispersion (normalized to ORIGINAL cs), matrix (each scale)
double *delta; // Vector for delta=ln(rho/rho0)-sigma^2/2
double *rho_old; //vector for previous values of rho
double *T_old; //vector for previous values of T
double *rho_T_field; //vector, conatining rho*T values, used for weighted averaging
complex *deltaK; //vector for Fourier transform of delta
complex *transform; //vector struct for Fourier transformation, see def. in FFT.h 
long *Xindex, *Yindex, *Zindex; //indices of coordinates
double *Xcoord, *Ycoord,*Zcoord; //coordinates of grid points
long *MinusKIndex; //vector, MinusKIndex[j] gives the index of -k[j]
double *Lvect; //scale sizes
long *NSpheres; //Number of spheres to use in Monte-Carlo for different scales
double *SphereVolume; //volumes of spheres of different sizes
double x,y,z; //relative coordinates of the MC spheres in pc
double xnorm,ynorm,znorm; //coordinates of the MC spheres in grid units
double *pointX,*pointY, *pointZ; //coordinates for subgrid
double SphereR, SphereMass, /*SphereRho,*/ SphereT; //parameters of MC spheres
double *rho_field_new, *T_field_new; // rho and T field in a collapsing subregions
double Vx_new, Vy_new, Vz_new; //relative velocity of subregion compared to ORIGINAL PARENT
double RelSphereR; //relative radius of subsphere to parent
double SphereMach;  //edge mah number of turbulent dispersion of subregion, using its own average c_s^2
double rhomax,Tmin, rholimit; // minimum of T used to calculate rholimit which is the minimum possible critical density at a scale
double MinSphereMass; //smallest mass that is enough to make the sphere self gravitate, used for MC sphere checks
// Stepsize control variables
double renormstepMax; // maximum of dt/tau among all modes
double dDeltaKMax; // maximum of sqrt(dt/tau*S) among all modes
long NfragmentsAlloc=100; //allocation size unit for cloud fragment number (unlikely to be used, as clouds don't have >100 fragments)
long Nfragments=0; //number of fragments this cloud has produced so far
double *FragmentMasses; //masses of substructures that fragmented out of this parent, used for virial calculation
double *FragmentDistances; //distances of substructures that fragmented out of this parent, used for virial calculation

//Control the efficency of heating
double heatingeff=1.0; //efficiency of heating, randomly determined, 100% default

//Shorthand variables
double Sigma0; //Original surface density
double *SqrtLvect; //square root of scales
double *InvLvectSq; //inverse square of scale sizes
double *Klength; //length of modes
double *InvKRelLength; //inverse relative k mode length, k_min/k
double *SqrtKRelLength; //square root of relative k mode length, sqrt(k/k_min)
double *OneOverKdensity4PiKcube; // shorthand for 1/(dNk/dVk(4*Pi*k^3))
double *renormstep; // shorthand for dt/tau where tau=sqrt(Lvect)
double *decaystep; //shorthand for 1-renormstep
double MachedgeSq; //edge Mach number square
double GammaMinusOne=gamma-1.0; //gamma-1

bool NotCMFSaved=true; //flag, tells us if this cloud is not yet included in the CMF file
bool ExistBeforeCMFTime=false; //flag, tells us if the cloud existed before CMF_Time. This means that it might still exist at CMF_Time

bool GammaLimitReached=false; //flag, tells if the cloud reached the limit in gamma when fragmentation is unlikely

bool HasFragmengted=false; //flag tells if the cloud has fragmented, used to check the reason for becoming not self gravitating

#if (SaveParentEvol==true)
char evolline[200]; //string to be saved, contains the global parameters of the parent cloud
char evolfilename[200]; //name of evolutionfile, which stores the time evolution of the global cloud
#endif

 /********************************************************************************************/
 	/* Initialize shorthands and variables */

#if (StochasticHeating==true)
if(StochasticHeating){
	//Set heating efficiency, lognormal(-S/2,S)~ exp(Gauss(0,1)*sqrt(S)-S/2), so the effective heating temperature is proportional to 0.25 power of this
	heatingeff=exp(0.25*(gaussrand()*sqrt(HeatingSVal)-HeatingSVal/2.0));
}
#endif

//Fix random number generator seed, if necessary for debugging
#if (RandomSeedFixed==true)
srand(0);
#endif

//  Get the individual process ID.
rank = MPI::COMM_WORLD.Get_rank();

#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"Initializing cloud variables... \n";
#endif

//Assign values to input variables
rho_field=inputstruct[0].rho_field;
T_field=inputstruct[0].T_field; //temperature field
//R0=inputstruct.R0; //initial size
X0=inputstruct[0].X0; Y0=inputstruct[0].Y0; Z0=inputstruct[0].Z0; //initial position
Machedge=inputstruct[0].Machedge;//initial Mach number
Vx=inputstruct[0].RelativeVelocityX; Vy=inputstruct[0].RelativeVelocityY; Vz=inputstruct[0].RelativeVelocityZ; //initial velocity
InitScale=inputstruct[0].InitScale; //relative scale below which turbulence is initialized, units of R0
tparent=inputstruct[0].tparent; //time at start

#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"Cloud variables from struct transferred, allocating memory for new variables... \n";
#endif


//Initialize modes
Lvect=vector_double(NScales);
LogScale(Lvect, 2.0*gridspacing, 0.85, NScales); //Minimum set to 2 grid spacing

//Scale variables
SqrtLvect=vector_double(NScales);
InvLvectSq=vector_double(NScales);
for ( i = 0; i < NScales; i++ ){
    SqrtLvect[i]=sqrt(Lvect[i]);
    InvLvectSq[i]=1.0/(Lvect[i]*Lvect[i]);
    //end for
    }
    
NSpheres=vector_long(NScales);//Number of spheres to use in Monte-Carlo for different scales
SphereVolume=vector_double(NScales);//volume of individual spheres
for ( i = 0; i < NScales; i++ ){
    NSpheres[i]=lround(NSphere0+MCFillFactor*pow(Lvect[i],-3.0));
    SphereVolume[i]=FourPiOver3*pow(Lvect[i],3.0); //Volumes in R units
    }

// Timesteps and shorthands
renormstep=vector_double(NgridTotal);
decaystep=vector_double(NgridTotal);
 
 //Allocating memory
 Mach=vector_double(NScales); //Mach number
 delta=vector_double(NgridTotal); //delta(x)=ln rho/rho0-variance^2
 S=vector_double(NgridTotal); //variance of delta(k)
 rho_old=vector_double(NgridTotal); //value for rho at previous timestep
 T_old=vector_double(NgridTotal); //value for T at previous timestep
 transform = (complex *) calloc( NgridTotal, sizeof(complex) ); //vector struct for Fourier transformation, see def. in FFT.h 
 deltaK=(complex *) calloc( NgridTotal, sizeof(complex) );  //delta(k), Fourier transform of delta
 rho_T_field=vector_double( NgridTotal);//vector, conatining rho*T values, used for weighted averaging
 pointX=vector_double(NgridTotal); //vector containing x coordinates for subgrid
 pointY=vector_double(NgridTotal); //vector containing y coordinates for subgrid
 pointZ=vector_double(NgridTotal); //vector containing z coordinates for subgrid
 rho_field_new=vector_double(NgridTotal); // density for subregion
 T_field_new=vector_double(NgridTotal); //T for subregion
 outputstruct=(struct inputstruct_s*) calloc(1,sizeof(struct inputstruct_s)); //allocate memory for output structure for recursion runs
 #if (DensitySave==true)
 	#if (ExtraVerboseOutput==true)
 		cout<<rank<<":\t"<<"Allocating density field struct...\n";
	#endif
	 fieldstruct=(struct densitystruct_s*) calloc(1,sizeof(struct densitystruct_s)); //allocate memory for field structure used to send data to master
 	#if (ExtraVerboseOutput==true)
 		cout<<rank<<":\t"<<"Finished allocating density field struct...\n";
	#endif
#endif
 FragmentMasses=vector_double(NfragmentsAlloc); // fragment masses, (it is unlikely for a cloud to produce more than a hundred
 FragmentDistances=vector_double(NfragmentsAlloc); // fragment masses, (it is unlikely for a cloud to produce more than a hundred

#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"Main allocation finsihed, setting up shorthands... \n";
#endif
 
 //Set up coordinate vectors
 Xindex=vector_long(NgridTotal); //x coordinate index
 Yindex=vector_long(NgridTotal); //y coordinate index
 Zindex=vector_long(NgridTotal); //z coordinate index
 Xcoord=vector_double(NgridTotal); //x coordinates on grid in R units
 Ycoord=vector_double(NgridTotal); //y coordinates on grid
 Zcoord=vector_double(NgridTotal); //z coordinates on grid
 
   #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG: Coordinate allocation done... \n";
#endif
 for ( i = 0; i < NgridTotal; i++ ){
     Xindex[i]= (i/Ngridsq); //x coordinate index
     Yindex[i]= ((i % Ngridsq)/Ngrid); //y coordinate index
     Zindex[i]= (i % Ngrid); //z coordinate index
     Xcoord[i]= Xindex[i]*gridspacing-1.0; //coordinate in R units with origin at center
     Ycoord[i]= Yindex[i]*gridspacing-1.0;
     Zcoord[i]= Zindex[i]*gridspacing-1.0;
     }//end for
  #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG: Index done... \n";
#endif
  // Set up vector containing indices of minus k vectors
  MinusKIndex=vector_long(NgridTotal);
  for ( i = 0; i < NgridTotal; i++ ){
  	j=(i/Ngridsq); //x coordinate
  	k=(i % Ngridsq)/Ngrid; //y coordinate
  	l=i % Ngrid; //z coordinate
  	//Flip them
  	j=(Ngrid-j)%Ngrid;
  	k=(Ngrid-k)%Ngrid;
  	l=(Ngrid-l)%Ngrid;
    MinusKIndex[i]=j*Ngridsq+k*Ngrid+l;
      }
      
  #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG: Minus index done... \n";
#endif
 
  //Set up k vector shorthands
 Klength=vector_double(NgridTotal); //length of k vectors
#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG:  Klength allocation  done... \n";
#endif
 InvKRelLength=vector_double(NgridTotal); //inverse relative k mode length, k_min/k
#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG:  InvKRelLengthallocation  done... \n";
#endif
 SqrtKRelLength=vector_double(NgridTotal); //relative square root k mode length, sqrt(k/k_min)
#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG:  SqrtKRelLength allocation  done... \n";
#endif
 OneOverKdensity4PiKcube=vector_double(NgridTotal); // shorthand for 1/(dNk/dVk(4*Pi*k^3))
#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG:  OneOverKdensity4PiKcube allocation  done... \n";
#endif
 for ( i = 0; i < NgridTotal; i++ ){
#if (ExtraVerboseOutput==true)
 if (MinusKIndex[i]>=NgridTotal){
 cout<<rank<<":\t"<<"DEBUG:  Too large MinusKindex for i="<<i<<" and minusindex="<<MinusKIndex[i]<<"... \n";
 }
#endif
 	if (i<=MinusKIndex[i]){ //each pair once time
 	
 		j=(i/Ngridsq); //x coordinate
  		k=(i % Ngridsq)/Ngrid; //y coordinate
  		l=i % Ngrid; //z coordinate
  		Klength[i]=sqrt((double)(j*j+k*k+l*l))*Kunit; //store length for this k vector
		//for -k index
		j=(MinusKIndex[i]/Ngridsq); //x coordinate
  		k=(MinusKIndex[i] % Ngridsq)/Ngrid; //y coordinate
  		l=MinusKIndex[i] % Ngrid; //z coordinate
  		Klength[MinusKIndex[i]]=sqrt((double)(j*j+k*k+l*l))*Kunit;
  	
  		//we only keep the shorter k, as large k corresponds to -k
  		if(Klength[i]<Klength[MinusKIndex[i]]){
  			Klength[MinusKIndex[i]]=Klength[i];
	  	}
	  	else{
	  		Klength[i]=Klength[MinusKIndex[i]];
	  	}
  		//calculate shorthands
    	InvKRelLength[i]=2.0*Kunit/Klength[i]; //inverse relative k mode length, k_(R)/k and Kunit=2*pi/(2 R) in R units
    	SqrtKRelLength[i]=1.0/sqrt(InvKRelLength[i]);
    	OneOverKdensity4PiKcube[i]=1.0/(4.0*M_PI*pow(Klength[i],3.0)*KSpaceDensity); // shorthand for 1/(dNk/dVk(4*Pi*k^3))
    	//for -k values, has to have the same values
    	InvKRelLength[MinusKIndex[i]]=InvKRelLength[i]; 
    	SqrtKRelLength[MinusKIndex[i]]=SqrtKRelLength[i];
    	OneOverKdensity4PiKcube[MinusKIndex[i]]=OneOverKdensity4PiKcube[i];
    }
 }
 
 /********************************************************************************************/
 	/* Initialize Starting values */
 	
#if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"Setting up initial values... \n";
#endif
     
 //Calculate total mass
 Mass=SumScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0); //get total mass in system
 MassInit=Mass; //initial mass
 MassWithinRadius=MassInit;
  #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG: Mass done... \n";
#endif

 //Original surface density
 Sigma0=Mass/(4.0*M_PI*R*R);

 //Calculate mean density
 RhoMean=AverageScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0); // it is more precise than using Mass/(geometric volume), due to discretization
 RhoMean0=RhoMean;
 #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG: Rhomean done... \n";
#endif
     
 //Initialize rho_T
 RhoTempFieldWithHeating( ( HeatingTemp (MassInit, R, heatingeff) ), T_field, rho_field, rho_T_field, NgridTotal); //calculate rho_T_field
  #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG: RhoT done... \n";
#endif
 //Init average temperatue by density averaging
 TMean=SumScalar_Sphere(rho_T_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0)/Mass;
  #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"DEBUG: Tmean done... \n";
#endif
 //Init gamma
 #if (IsothermalEOS==true)
 gamma=1.0;
 #else
 	 gamma=gammafunc(RhoMean0, TMean);
 #endif
//Enforce virial equilibrium
MachedgeSq=(Gconst*MassWithinRadius/R)/SoundSpeedSquare(TMean, gamma)*QVirial-1.0;
if(MachedgeSq<0){ //The cloud is thermally supported!
	MachedgeSq=0.0;
}
 #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"Virial equilibrium enforced, T="<<TMean<<". Original Mach="<<Machedge<<"\t New Mach="<<sqrt(MachedgeSq)<<"\n";
 #endif

 
 #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"Starting calculation with M="<<Mass<<", T="<<TMean<<", R="<<R<<", Mach="<<sqrt(MachedgeSq)<<", Vx="<<Vx<<", Vy="<<Vy<<", Vz="<<Vz<<"\n";
 #endif

//Initialize deltaK
    for ( j = 0; j < NgridTotal; j++ ){
    	rho_field[j]=rho_field[j];
        deltaK[j].Re=log(rho_field[j]/RhoMean0); //real part, this is just delta(x) at t=0
                                      }
    //Fourier transform of density field
    fft3D(deltaK, Ngrid, Ngrid, Ngrid, +1); // Fourier transform of delta(x), NO NORMALIZATION!
    //Applying normlaization, multiplying by dxdydz=(2R/N)^3
    for ( j = 0; j < NgridTotal; j++ ){
        deltaK[j].Re=deltaK[j].Re/NgridTotal; 
        deltaK[j].Im=deltaK[j].Im/NgridTotal; 
                                      }//endfor
   
  //Initialize variables used for cloud contraction calculation
  y0=vector_double(dimension); //this will have the results
  y1=vector_double(dimension); //these are used to check for error
  y2=vector_double(dimension);
  y0[0]=1.0; //y[0] - relative radius
  y0[1]=Sigma0;//y[1] - surface density, sigma=M/(4 Pi R^2) 
  //MachedgeSq=Machedge*Machedge; //will be used for y[2] - current edge Mach number squared
  y0[2]=MachedgeSq; //Mach squared
  y0[3]=1.0; //initial mass ratio, M/M_init
  
  //allocating memory for results
  RelSizeVector=vector_double(MaxStepNum); //relative size
  RelSizeVector[0]=y0[0];
  //Setting alternates
  y1[0]=y0[0];y1[1]=y0[1];y1[2]=y0[2];y1[3]=y0[3];
  y2[0]=y0[0];y2[1]=y0[1];y2[2]=y0[2];y2[3]=y0[3];


/********************************************************************************************/
	/*  Turbulence Initialization */
	
// GMCs should start out having fully developed turbulence, which means that the homogeneous initial conditions are not good.
// This is a copy of the collapse time evolution where no collapse occurs, only the rho and T fields evolve in time
// Since the parent can only initialize turbulence up to its resolution, it is initialized here for smaller scales

if((InitScale>0) && (Mass>=MinMass) && (MachedgeSq>0.0)){ //there is a small enough scale that is not initialized by the parent routine and we want to actually evolve this
	t=0;
	//Gamma polytorpic index as time evolves
	#if(IsothermalEOS)
	gamma=1.0;	
	#else
		gamma=gammafunc(RhoMean0, TMean);
	#endif
	GammaMinusOne=gamma-1.0;

    // Caculate Mach number, according using turbulent scaling Mach^2~R
    for ( j = 0; j < NScales; j++ ){
        Mach[j]=SqrtLvect[j]*Machedge;
    }//end for

	//crossing time in Myr
    tcrossing=2*R/sqrt((QVirial*Gconst*MassInit/R0)-SoundSpeedSquare(TMean, gamma))*PcPerMyr;
    timeconversionfactor=2.0*pow(QVirial,-1.5)*pow(R0,1.5)/sqrt(MassInit)*TimeUnit; //freefall time of the cloud, in Myr
	timescalecorrection=timeconversionfactor/tcrossing;
    
    // Theoretical maximum of stepping among different modes (max of dt/tau and sqrt(S*dt))
    renormstepMax=dt*timescalecorrection*sqrt((double )Ngrid/2); // dt/tau=dt/t0*sqrt(lambda/R) so we get max for highest k mode, kmax=2*R/Ngrid
    dDeltaKMax=sqrt(dt*timescalecorrection*log(1.0+2*BSq*MachedgeSq)*0.0562698);  // S*dt/tau=dt*log(1+2*0.25*Machedge^2)/(4.0*M_Pi*sqrt(2))) for largest mode (lambda=2*R)
    //Check if dt is too big or too small and change accordingly
    dt=SetOptimaldt(dt, renormstepMax, dDeltaKMax);

    //Calculate S(k)
    for ( j = 1; j < (NgridTotal); j++ ){ //we start from 1, because j=0 is the K=0 constant mode
        renormstep[j]=dt*timescalecorrection*SqrtKRelLength[j];   //dt/tau where tau is proportional to sqrt(Lvect)
        decaystep[j]=1.0-renormstep[j];     //shorthand for 1-renormstep
		S[j]=log(1.0+BSq*MachedgeSq*InvKRelLength[j])*OneOverKdensity4PiKcube[j]; //S=log(1+b^2*Mach^2)/( 4 Pi k^3 dN/dV) where we used Mach^2~R
        }
    
	while(t<InitTime){
    	NormFactor=0.0;
    	//Calculate the variance of delta(k) and evolve deltaK
        //go only till half of the X values, since we will get the rest from delta(k)*=delta(-k)
    	for ( j = 1; j < (NgridTotal); j++ ){ //we start from 1, because j=0 is the K=0 constant mode
    	  if(InvKRelLength[j]<InitScale){ //means that this scale is not initialized by the parent
			//Evolve deltaK
			if (j<MinusKIndex[j]){ //different indices, first time
            	//Real part
            	temp=gaussrand()*sqrt(S[j]*renormstep[j]);//gaussian random*sqrt(2 dS dt/tau),
            	deltaK[j].Re= deltaK[j].Re*decaystep[j]+temp; //Fokker-Planck evolution following Eq. 54, assuming S_re=1/2 S
        		//Set up delta(-k) to be conjugate, so that final result is real, this assumes Ngrid is even
        		deltaK[MinusKIndex[j]].Re=deltaK[MinusKIndex[j]].Re*decaystep[j]+temp;
            	//Imaginary part
        		temp=gaussrand()*sqrt(S[j]*renormstep[j]);//gaussian random*sqrt(2 dS dt/tau),
        		deltaK[j].Im= deltaK[j].Im*decaystep[j]+temp; //Fokker-Planck evolution following Eq. 54, assuming S_im=1/2 S
        		//Set up delta(-k) to be conjugate, so that final result is real, this assumes Ngrid is even
        		deltaK[MinusKIndex[j]].Im=deltaK[MinusKIndex[j]].Im*decaystep[j]-temp;
        	}
        	else{ 
        		if (j==MinusKIndex[j]){ //it does not have imaginary values, this is the only way it is equal to its conjugate
            		temp=gaussrand()*sqrt(S[j]*renormstep[j]);//uniform random [0;1]*(variance),
            		deltaK[j].Re= deltaK[j].Re*decaystep[j]+temp; //evolution following Eq. 54, assuming S_re=1/2 S
        		}
        	}
		  }//end if (init)
        } //end for
    
    	// Calculate real space delta field
          //Prepare data
    	for ( j = 0; j < NgridTotal; j++ ){
        	transform[j].Re=deltaK[j].Re; //real part
        	transform[j].Im=deltaK[j].Im; //Imaginary part
        }//endfor

          //Fourier transform of deltaK
    	fft3D(transform, Ngrid, Ngrid, Ngrid, -1); // Inverse Fourier transform (back to real space), normalization applied at next step
    
          //Get real space log density
    	for ( j = 0; j < NgridTotal; j++ ){
        	delta[j]=transform[j].Re; //real part only, no imaginary density
        	//Get real space density
        	rho_old[j]=rho_field[j]; //storing previous result
        	rho_field[j]=exp(delta[j]); // gets rho without normalization
        }//endfor
    	
      // Normalize rho; ln(rho/rho0)=delta-sigma^2/2 => rho=exp(delta)*NormFactor
        //NormFactor= Mass/(Integral of non normalized rho)
    	NormFactor=Mass/SumScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0);
    	for ( j = 0; j < NgridTotal; j++ ){
        	rho_field[j]*=NormFactor; // normalizes rho to have the correct mass in the collapsing sphere
        	if (rho_field[j]<RhoEmpty){ //should have at least this value
        		rho_field[j]=RhoEmpty;
			}
    	}
    
    	//Evolve temperature field
    	#if(IsothermalEOS==false)
	    for ( j = 0; j < NgridTotal; j++ ){
	        T_old[j]=T_field[j];
	        //Evolve T locally as a polytrope
	    	T_field[j]=T_field[j]*pow(rho_field[j]/rho_old[j],GammaMinusOne);
	        //rho_T_field[j]=T_field[j]*rho_field[j];
	    }
        #endif

    	RhoTempFieldWithHeating( ( HeatingTemp (MassInit, R, heatingeff) ), T_field, rho_field, rho_T_field, NgridTotal); //calculate rho_T_field
		t+=dt;	
	}
    //Get mean temperature on largest scale, weighted average
    TMean=SumScalar_Sphere(rho_T_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0)/Mass;
}

//Initialize gamma for the isothermal case
#if (IsothermalEOS)
	gamma=1.0;
	GammaMinusOne=0.0;
#endif

//Enforce virial equilibrium after turbulence initialization (EOS and other stuff might mess with it)
MachedgeSq=(Gconst*MassWithinRadius/R)/SoundSpeedSquare(TMean, gamma)*QVirial-1.0;
if(MachedgeSq<0){ //The cloud is thermally supported!
	MachedgeSq=0.0;
}
 #if (ExtraVerboseOutput==true)
 cout<<rank<<":\t"<<"Virial equilibrium enforced, T="<<TMean<<". Original Mach="<<Machedge<<"\t New Mach="<<sqrt(MachedgeSq)<<"\n";
 #endif

                           

/************************************************************************************************************/
   /*  Time evolution (Collapse)  */ 
i=0; //reset i
t=0.0;//reset t
timeconversionfactor=2.0*pow(QVirial,-1.5)*pow(R0,1.5)/sqrt(MassInit)*TimeUnit; //freefall time of the cloud, in Myr
CMF_Time=FindCMFTime(tcurrent); //find closest observation time
CMF_Time_Index=FindCMFTimeIndex(tcurrent); //find closest observation time

//Check current time against CMF observation time
if(CMF_Time_Index!=-1 && tcurrent<=CMF_Time){
	ExistBeforeCMFTime=true; //this cloud existed before the observed time, question is: does it survive till that?
}

//repeat until we get to small enough scale, or we run out of steps, or too small fragment to matter (at the start) or no more turbulence to drive anything
while ((y0[0]>ymin) && (t<Tmax) && (Mass>=MinMass || t>0) && (MachedgeSq>0.0)){ 

if (i<FIRSTTIMESTEPNUM){ //the first few time steps are very small so that the hierarchy for newly initialized scales can be established
	dt=FIRSTTIMESTEP;
}
else{
	dt=DEFAULTTIMESTEP; //use default time steps and reduce it if necessary
	// Theoretical maximum of stepping among different modes (max of dt/tau and sqrt(S*dt))
    renormstepMax=dt*timescalecorrection*sqrt((double )Ngrid/2); // dt/tau=dt/t0*sqrt(lambda/R) so we get max for highest k mode, kmax=2*R/Ngrid
    dDeltaKMax=sqrt(dt*timescalecorrection*log(1.0+2*BSq*MachedgeSq)*0.0562698);  // S*dt/tau=dt*log(1+2*0.25*Machedge^2)/(4.0*M_Pi*sqrt(2))) for largest mode (lambda=2*R)
    //Check if dt is too big or too small and change accordingly
    dt=SetOptimaldt(dt, renormstepMax, dDeltaKMax);
}

//Find the appropriate timestep, using RK4
  do{ 
     //Setting alternates (y(2) does not change)
     y0[2]=MachedgeSq; //update edge Mach number squared (recalculated at the end of rho and T field evolution step)
     y0[3]=MassWithinRadius/MassInit;
	 y1[0]=y0[0];
     y1[1]=y0[1];
     y1[2]=y0[2];
     y1[3]=y0[3];
     y2[0]=y0[0];
     y2[1]=y0[1];
     y2[2]=y0[2];
     y2[3]=y0[3];
    // Take a step with dt   
    rk4_step ( y1, t, dt, dimension);
    // Take two steps with dt/2
    rk4_step ( y2, t, dt/2.0, dimension);
    rk4_step ( y2, t, dt/2.0, dimension);
    //Compare
    err=fabs(y1[0]-y2[0]); //absolute error
    if ( (dt>MinStep) && ((err>MaxAbsErr)  || (y2[0]<=0.0) ) ){ //reduce step size if error is too big or negative value
       dt=dt/2.0; //decreases step size
       }
    //end do-while loop
    }while( (dt>MinStep) && ((err>MaxAbsErr)  || (y2[0]<=0.0) ) ); //reduce step size if error is too big or negative value
    //Take step
    ++i; //increases step number
	y0[0]=y2[0]; //take value
	//Store results
    RelSizeVector[i]=y0[0]; //relative size at the time
    //Get new sigma
    y0[1]=Sigma0/(RelSizeVector[i]*RelSizeVector[i]);
    
//Time to calculate parameters at this new size

    //Size at time t
    R=R0*RelSizeVector[i];
    
    //Calculate the amount of mass enclosed in the cloud (inlcuding fragments)
    MassWithinRadius=Mass;
    for(j=0;j<Nfragments;j++){
    	if(FragmentDistances[j]<R){
    		MassWithinRadius+=FragmentMasses[j];
		}
	}

    //Mean density at this time
    RhoMean=RhoMean0/pow(RelSizeVector[i],3.0);

    //Gamma polytorpic index as time evolves
    #if (IsothermalEOS==false)
    	gamma=gammafunc( RhoMean, TMean );
    	GammaMinusOne=gamma-1.0;
	#endif
	
	//crossing time in Myr
    tcrossing=2*R/sqrt((QVirial*Gconst*MassWithinRadius/R)-SoundSpeedSquare(TMean, gamma))*PcPerMyr;
	timescalecorrection=timeconversionfactor/tcrossing;
    
    // Caculate Mach number, according using turbulent scaling Mach^2~R for subregions of different sizes
    for ( j = 0; j < NScales; j++ ){
        Mach[j]=SqrtLvect[j]*sqrt( MachedgeSq );
    }//end for
    NormFactor=0.0; //reset
    //Calculate S(k), the variance of delta(k) and evolve deltaK
          //go only till half of the X values, since we will get the rest from delta(k)*=delta(-k)
    for ( j = 1; j < (NgridTotal); j++ ){ //we start from 1, because j=0 is the K=0 constant mode
        renormstep[j]=dt*SqrtKRelLength[j]*timescalecorrection;   //dt/tau where tau=crossing time=sqrt(Lvect)*(crossing time on cloud scale)
		decaystep[j]=1.0-renormstep[j];     //shorthand for 1-renormstep
		S[j]=log(1.0+BSq*MachedgeSq*InvKRelLength[j])*OneOverKdensity4PiKcube[j]; //S=log(1+b^2*Mach^2)/( 4 Pi k^3 dN/dV) where we used Mach^2~R
		//Evolve deltaK
		if (j<MinusKIndex[j]){ //different indices, first time
            //Real part
            temp=gaussrand()*sqrt(S[j]*renormstep[j]);//gaussian random*sqrt(2 dS dt/tau),
            deltaK[j].Re= deltaK[j].Re*decaystep[j]+temp; //evolution following Eq. 54, assuming S_re=1/2 S
        	//Set up delta(-k) to be conjugate, so that final result is real, this assumes Ngrid is even
        	deltaK[MinusKIndex[j]].Re=deltaK[MinusKIndex[j]].Re*decaystep[j]+temp;
            //Imaginary part
        	temp=gaussrand()*sqrt(S[j]*renormstep[j]);//gaussian random*sqrt(2 dS dt/tau),
        	deltaK[j].Im= deltaK[j].Im*decaystep[j]+temp; //evolution following Eq. 54, assuming S_im=1/2 S
        	//Set up delta(-k) to be conjugate, so that final result is real, this assumes Ngrid is even
        	deltaK[MinusKIndex[j]].Im=deltaK[MinusKIndex[j]].Im*decaystep[j]-temp;
        }
        else{ 
        	if (j==MinusKIndex[j]){ //it does not have imaginary values, this is the only way it is equal to its conjugate
            	temp=gaussrand()*sqrt(S[j]*renormstep[j]);//uniform random [0;1]*(variance),
            	deltaK[j].Re= deltaK[j].Re*decaystep[j]+temp; //evolution following Eq. 54, assuming S_re=1/2 S
        	}
        }
 
        } //end for
    
    // Calculate real space delta field
          //Prepare data
    for ( j = 0; j < NgridTotal; j++ ){
        transform[j].Re=deltaK[j].Re; //real part
        transform[j].Im=deltaK[j].Im; //Imaginary part
                                      }//endfor

          //Fourier transform of deltaK
    fft3D(transform, Ngrid, Ngrid, Ngrid, -1); // Inverse Fourier transform (back to real space), normalization applied at next step
    
          //Get real space log density
    for ( j = 0; j < NgridTotal; j++ ){
        delta[j]=transform[j].Re; //real part only, no imaginary density
    	//Get real space density
        rho_old[j]=rho_field[j]; //storing previous result
        rho_field[j]=exp(delta[j]); // gets rho without normalization
        }
        
        rhomax=RhoEmpty; //set it to default value, will find the actual value later
      // Normalize rho; ln(rho/rho0)=delta-sigma^2/2 => rho=exp(delta)*NormFactor
        //NormFactor= Mass/(Integral of non normalized rho)
    NormFactor=Mass/SumScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0);
    for ( j = 0; j < NgridTotal; j++ ){
        rho_field[j]*=NormFactor; // normalizes rho to have the correct mass in the collapsing sphere
        if (rho_field[j]<RhoEmpty){ //should have at least this value
        	rho_field[j]=RhoEmpty;
		}
		//Find maximum of rho
		if (rho_field[j]>rhomax)
		{
			rhomax=rho_field[j]; //sets the new maximum
		}
    }
    //IMPORTANT: This rho includes the increase due to contraction!!!!!!
    
    //Evolve temperature field
    Tmin=TMean; //set it as a default value
    #if(IsothermalEOS==false)
	for ( j = 0; j < NgridTotal; j++ ){
	    T_old[j]=T_field[j];
	    //Evolve T locally as a polytrope
	    T_field[j]=T_field[j]*pow(rho_field[j]/rho_old[j],GammaMinusOne);
	    //rho_T_field[j]=T_field[j]*rho_field[j];
	    if(Tmin>T_field[j]){
	        Tmin=T_field[j]; //find smallest T
		}
	}
	#endif
	//calculate rho_T_field
    RhoTempFieldWithHeating( ( HeatingTemp (MassWithinRadius, R, heatingeff) ), T_field, rho_field, rho_T_field, NgridTotal); 
    //Get mean temperature on largest scale, weighted average
    TMean=SumScalar_Sphere(rho_T_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0)/Mass;
    
    //Evolve edge Mach number so that virial equilibrium is kept, use initial mass, as the fragments are still present and contributing
    MachedgeSq=(Gconst*MassWithinRadius/R)/SoundSpeedSquare(TMean, gamma)*QVirial-1.0;
    if(MachedgeSq<0){ //The cloud is thermally supported!
    	MachedgeSq=0.0;
	}
	
	//Save eolution of global parameters of the parent cloud
#if (SaveParentEvol==true)
	if(X0==0 && Y0==0.0 && Z0==0.0){
		//Saved variables: renorm time, real time, radius, relative size, mass, mean temperature, gamma, Mach number, sound speed, sonic length, sonic mass
		sprintf( evolline, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", t,tcurrent,R,RelSizeVector[i],Mass,TMean,gamma,sqrt(MachedgeSq),sqrt(SoundSpeedSquare(TMean, gamma)),R/MachedgeSq, (2.0/QVirial*SoundSpeedSquare(TMean, gamma)*R/MachedgeSq/Gconst ) );
		sprintf( evolfilename, "M_%g_globalevol.txt", MassInit );
		WriteLineToTextFile(evolline, evolfilename);
	}

#endif
  
    //Check if there is a collapsing subregion by throwing down spheres and checking if they are self-gravitating
#if (StartWithLargeFragments==true)
	for ( j=(NScales-1); j>=0; j-- ){ //for each scale, start at the largest (this creates a bit of an overhead but the small scales are initialized in a more consistent manner)
#else
	for ( j=0; j<NScales; j++ ){ //for each scale, start looking for fragments at the smallest scales
#endif
        SphereR=R*Lvect[j];
        //minimum possible rho critical for this scale
        rholimit=RhoCrit(Tmin,SphereR, (Mach[j]*Mach[j]), gamma);
        //Minimum mass to have self gravitating cloud of this size
        MinSphereMass=rholimit*SphereVolume[j]*pow(R,3.0);
      if(FragmentAllowed && rhomax>rholimit){ //with max rho and min T it could collapse
	  for ( k = 0; k < NSpheres[j]; k++ ){ //each Monte-Carlo sphere
      //Get center of sphere (random)
      do{
       	x=((double)rand()/((double)RAND_MAX + 1.0)-0.5)*2.0*(R-SphereR); //x coordinate between -R and R
       	y=((double)rand()/((double)RAND_MAX + 1.0)-0.5)*2.0*(R-SphereR); //y coordinate
      	z=((double)rand()/((double)RAND_MAX + 1.0)-0.5)*2.0*(R-SphereR); //z coordinate
      }while( x*x+y*y+z*z> (R-SphereR)*(R-SphereR) ); //repeat until we get a sphere which is inside
      xnorm=x/R;ynorm=y/R;znorm=z/R; //coordinates in grid units
      //Get mean density in sphere by taking total mass in the sphere divided by volume
      SphereMass=SumScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R,SphereR , x, y, z); //mass in sphere
      if(SphereMass>MinSphereMass && SphereMass>IgnoreMass && SphereMass<(MaxFragmentRelativeMass*Mass) ){ //at least the required mass is in it and large enough that it is not ignored but not too large comaped to parent
        // Get mean temperature in sphere
        SphereT=SumScalar_Sphere(rho_T_field, Ngrid,(-R), R, (-R), R, (-R), R, SphereR, x, y, z)/SphereMass;
        //Find the edge Mach number normalized to its own temperature, using the fact that Mach~1/c_s~1/sqrt(T)
        SphereMach=Mach[j]*sqrt(TMean/SphereT);

        //using Jeans condition to see if it is self-gravitating
        if( IsSelfGravitating(SphereMass, SphereT, SphereR, (SphereMach*SphereMach), gamma, StrictRhoFactor) )
        {
		//Check if it will continue to collapse after we zoom in, so let's give a higher estimate for the Mach number, skipped by default
		if(SkipExtraCollapseCheck==true || 0<((SphereMass/SphereR*Gconst)/SoundSpeedSquare(HeatingTemp (SphereMass, SphereR, 1.0)+0.1, gamma)*QVirial-1.0)){ //the 0.1 is just to make sure we never get zero (e.g. case of no heating)
		 	//Tell the Master that the we have data for a new run and a recurison is needed
			message=RecursionMessage; //message= "Need to call a recursion. Sending input soon"
			mpierror=MPI_Ssend( &message, 1, MPI_INT, 0, SendMessageTag, MPI_COMM_WORLD);//send message
			
            //Interpolate fields
            for ( l = 0; l < NgridTotal; l++ ){ 
                //Coordinate setup
                pointX[l]=x+SphereR*Xcoord[l]; //x coordinate for point l in subregion
                pointY[l]=y+SphereR*Ycoord[l]; //y coordinate for point l in subregion
                pointZ[l]=z+SphereR*Zcoord[l]; //z coordinate for point l in subregion
                }
            //Set new rho to empty values, we will fill out the rest afterwards
            vector_set( rho_field_new, RhoEmpty, NgridTotal );
            //Go again for interpolation
            for ( l = 0; l < NgridTotal; l++ ){
                //Interpolate values for points in sphere
                rho_field_new[l]=trilinear_interp(pointX[l],pointY[l],pointZ[l],rho_field,Ngrid,Ngrid,Ngrid,(-R), R, (-R), R, (-R), R,0,RhoEmpty);
                T_field_new[l]=trilinear_interp(pointX[l],pointY[l],pointZ[l],T_field,Ngrid,Ngrid,Ngrid,(-R), R, (-R), R, (-R), R,0,TEmpty);
                }
            //Renormalize subgrid to fix errors introduced by finite grid, this way mass is conserved
            NormFactor=SphereMass/SumScalar_Sphere(rho_field_new, Ngrid,(-SphereR), SphereR, (-SphereR), SphereR, (-SphereR), SphereR,SphereR , 0.0, 0.0, 0.0);
            vector_const_multiply( rho_field_new, rho_field_new, NormFactor, NgridTotal );
            //Remove density from original grid
            RelSphereR=Lvect[j];//using R=1 units, so RelSphereR=Lvect[j]
            for ( l = 0; l < NgridTotal; l++ ){
            //we will go through all point in the parent sphere but we only modify those within the collapsing sphere, 
                 if ( (Xcoord[l]-xnorm)*(Xcoord[l]-xnorm) + (Ycoord[l]-ynorm)*(Ycoord[l]-ynorm) + (Zcoord[l]-znorm)*(Zcoord[l]-znorm) <= RelSphereR*RelSphereR )
				 {
					rho_field[l]=RhoEmpty;// set density as empty
                    T_field[l]=TEmpty;// set tempreature as the mean
                 }//end if
                }//end for loop through all grid points
			//Renormalize parent to fix errors inttroduced by finite grid, this way mass is conserved
			NormFactor=(Mass-SphereMass)/SumScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0);
			vector_const_multiply( rho_field, rho_field, NormFactor, NgridTotal );
                
            //Velocity for subgrid relative to ORIGINAL PARENT; velocities add up
            temp=sqrt(MachedgeSq/3.0*sqrt((x*x+y*y+z*z)/(R*R))*SoundSpeedSquare(TMean, gamma)); //velocity dispersion
            #if (ExtraVerboseOutput==true)
            cout<<rank<<":\t"<<"Velocity dispersion (1D) at substructure distance: "<<temp<<" m/s\n";            
            #endif
            Vx_new=Vx+gaussrand()*temp;
            Vy_new=Vy+gaussrand()*temp;
            Vz_new=Vz+gaussrand()*temp;
			        
            // We have the subgrid, time to call the recursion
            
            HasFragmengted=true; //the cloud has fragmented
            
            FragmentMasses[Nfragments]=SphereMass; //mass of substructure, in Msolar
			FragmentDistances[Nfragments]=sqrt(x*x+y*y+z*z); //distances of substructure, in pc
			Nfragments++;
			//Check if we have don't have enogh allocations and realloc
			if(NfragmentsAlloc==Nfragments){
				FragmentMasses=(double*) realloc(FragmentMasses, (NfragmentsAlloc+100)* sizeof(double));
				FragmentDistances=(double*) realloc(FragmentDistances, (NfragmentsAlloc+100)* sizeof(double));
				NfragmentsAlloc=NfragmentsAlloc+100;
				if(FragmentMasses==NULL || FragmentDistances==NULL){ 
					cout<<"Memory allocation failure for FragmentMasses! Exiting!";
					exit(EXIT_FAILURE);
				}
			}
			
			//Transfer to output structure
			for ( l = 0; l < NgridTotal; l++ ){
				outputstruct[0].rho_field[l]=rho_field_new[l];
				outputstruct[0].T_field[l]=T_field_new[l];
			}
			outputstruct[0].InitScale=InitAllowed*gridspacing/Lvect[j]; //smallest relative scale of the subregion not initialized
			outputstruct[0].Machedge=SphereMach;
			outputstruct[0].R0=SphereR;
			outputstruct[0].X0=x+X0;outputstruct[0].Y0=y+Y0;outputstruct[0].Z0=z+Z0;//we need to be careful with the coordinates, the absolute coordinate=relative to direct parent + coordinate of parent
			outputstruct[0].RelativeVelocityX=Vx_new;outputstruct[0].RelativeVelocityY=Vy_new;outputstruct[0].RelativeVelocityZ=Vz_new;
			outputstruct[0].tparent=tcurrent;
			
			#if (ExtraVerboseOutput==true)
			cout<<rank<<":\t"<<"Evolvefield found new structure with parameters: M="<<SphereMass<<", T="<<SphereT<<", R="<<SphereR<<", Mach="<<SphereMach<<", Vx="<<Vx_new<<", Vy="<<Vy_new<<", Vz="<<Vz_new<<"\n";
			#endif
			
			//Send output structure to Master
			mpierror=MPI_Ssend( outputstruct, 1, mpi_inputstruct_type, 0, SendInputTag, MPI_COMM_WORLD);
			
			#if (VerboseOutput==true)
			cout<<rank<<":\t"<<"Evolefield sent data\n";
			#endif
  
		//Now we have a new T and rho fields, time to update the parameters
            //Get total mass
            Mass=SumScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0);
            //Less mass, which changes Sigma0
            Sigma0=Mass/(4.0*M_PI*R0*R0);
            //New gamma
            #if (IsothermalEOS==false)
            gamma=gammafunc(RhoMean, TMean );
            #endif
            //We need to update RhoMean too
            RhoMean=AverageScalar_Sphere(rho_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0);
            RhoMean0=RhoMean*pow(RelSizeVector[i],3.0); //new value of "initial" mean density, so that the scaling of rhomean=rhomean0*size^(-3) still works
            //Set up rho_T field
            RhoTempFieldWithHeating( ( HeatingTemp (MassWithinRadius, R, heatingeff) ), T_field, rho_field, rho_T_field, NgridTotal); //calculate rho_T_field
            //Get mean temperature on largest scale, weighted average
    		TMean=SumScalar_Sphere(rho_T_field, Ngrid,(-R), R, (-R), R, (-R), R, R, 0.0, 0.0, 0.0)/Mass;
    		//Evolve edge Mach number so that virial equilibrium is kept
    		MachedgeSq=(MassWithinRadius/R*Gconst)/SoundSpeedSquare(TMean, gamma)*QVirial-1.0;
    		if(MachedgeSq<0){ //The cloud is thermally supported!
    			MachedgeSq=0.0;
			}
            //We have lost mass, that will influence the collapse of the remaining material, so we update sigma for the RK4 solver
            y0[1]=Mass/(4.0*M_PI*R*R);y1[1]=y0[1];y2[1]=y0[1];
            //We need to reset deltaK, using same method as during initialization
               //Filling it up with log rho values
            for ( l = 0; l < NgridTotal; l++ ){
                deltaK[l].Re=log(rho_field[l]/RhoMean); //real part, this is just delta(x) at t=0
                deltaK[l].Im=0.0;  //Log density is real
                                      }
            //Fourier transform of density field
            fft3D(deltaK, Ngrid, Ngrid, Ngrid, +1); // Fourier transform of delta(x) to get delta(k)
            //Normalization of modes, as the FFT has no normalization
            for ( l = 0; l < NgridTotal; l++ ){
        		deltaK[l].Re=deltaK[l].Re/NgridTotal; 
        		deltaK[l].Im=deltaK[l].Im/NgridTotal; 
                }//endfor

        //Now we have implemented the changes to rho, T and are ready to continue the time evolution

            }//end collapse if (Mach number check)
        	}//end fragment collapse check if
		  } //end if for spheremass
        }// end for cycle for MC spheres
    }//end if for sphere size
    }// end for cycle for MC sphere sizes

  
  /*******************************/
  // CMF evolution saving
  if(CMFEvolSave && ExistBeforeCMFTime && NotCMFSaved && tcurrent>=CMF_Time){ //exists at the desired time and we haven't saved it yet
    #if (VerboseOutput==true)
  	cout<<rank<<":\t"<<"Slave found a cloud\n";
  	#endif
	//Tell the Master it found a core and transfer the data
	message=CloudFoundMessage; //message= "I found a self-gravitating cloud"
	mpierror=MPI_Ssend( &message, 1, MPI_INT, 0, SendMessageTag, MPI_COMM_WORLD);//send message
	//Start sending data
		//Send time index
	mpierror=MPI_Ssend( &CMF_Time_Index, 1, MPI_LONG, 0, SendTimeIndexTag, MPI_COMM_WORLD);
		//cloud parameters
	mpierror=MPI_Ssend( &Mass, 1, MPI_DOUBLE, 0, SendMassTag, MPI_COMM_WORLD); //mass
	mpierror=MPI_Ssend( &X0, 1, MPI_DOUBLE, 0, SendXTag, MPI_COMM_WORLD); //X
	mpierror=MPI_Ssend( &Y0, 1, MPI_DOUBLE, 0, SendYTag, MPI_COMM_WORLD); //Y
	mpierror=MPI_Ssend( &Z0, 1, MPI_DOUBLE, 0, SendZTag, MPI_COMM_WORLD); //Z
	mpierror=MPI_Ssend( &Vx, 1, MPI_DOUBLE, 0, SendVXTag, MPI_COMM_WORLD); //VX
	mpierror=MPI_Ssend( &Vy, 1, MPI_DOUBLE, 0, SendVYTag, MPI_COMM_WORLD); //VY
	mpierror=MPI_Ssend( &Vz, 1, MPI_DOUBLE, 0, SendVZTag, MPI_COMM_WORLD); //VZ
		//extra parameters
	mpierror=MPI_Ssend( &MassInit, 1, MPI_DOUBLE, 0, SendMassInitTag, MPI_COMM_WORLD); //initial mass
	temp=sqrt(MachedgeSq);
	mpierror=MPI_Ssend( &temp, 1, MPI_DOUBLE, 0, SendMachTag, MPI_COMM_WORLD); //mach number
	mpierror=MPI_Ssend( &TMean, 1, MPI_DOUBLE, 0, SendTemperatureTag, MPI_COMM_WORLD); //Temperature
	mpierror=MPI_Ssend( &RhoMean, 1, MPI_DOUBLE, 0, SendRhoTag, MPI_COMM_WORLD); //average density
	mpierror=MPI_Ssend( &R, 1, MPI_DOUBLE, 0, SendRadiusTag, MPI_COMM_WORLD); //radius
	
	#if (DensitySave==true) //send fields to master who will save them into a file  
	//Tell the Master it found a the fields at the right time
	message=FieldsMessage; //message= "I found a self-gravitating cloud"
	mpierror=MPI_Ssend( &message, 1, MPI_INT, 0, SendMessageTag, MPI_COMM_WORLD);//send message
	//Send time index
	mpierror=MPI_Ssend( &CMF_Time_Index, 1, MPI_LONG, 0, SendTimeIndexTag, MPI_COMM_WORLD);
		//Transfer data to field structure
	for ( l = 0; l < NgridTotal; l++ ){
		fieldstruct[0].rho_field[l]=rho_field[l];
		fieldstruct[0].T_field[l]=T_field[l];
		fieldstruct[0].xcoord[l]=X0+R*Xcoord[l];
		fieldstruct[0].ycoord[l]=Y0+R*Ycoord[l];
		fieldstruct[0].zcoord[l]=Z0+R*Zcoord[l];
	}
	fieldstruct[0].volume=FourPiOver3*pow(R,3.0)/((double)NgridTotal); //volume of a single cell
		//Send density structure to Master
		mpierror=MPI_Ssend( fieldstruct, 1, mpi_densitystruct_type, 0, SendFieldsTag, MPI_COMM_WORLD);
	#endif

	#if (VerboseOutput==true)
	cout<<rank<<":\t"<<"Slave sent a cloud\n";
	#endif
	
    NotCMFSaved=false; //mark that it is saved
  }
  
  //Update CMF time
	if(CMF_Time!=-1.0 && CMF_Time<tcurrent){ //if there is a valid CMF time check if we have gone past it
    	CMF_Time=FindCMFTime(tcurrent); //find next closest observation time
    	CMF_Time_Index=FindCMFTimeIndex(tcurrent); //find closest observation time
    	if(CMF_Time==-1.0){
    		ExistBeforeCMFTime=false; //no valid CMFTime target, flag is meaningless
		}
		else{
			ExistBeforeCMFTime=true; //exists before observation time
			NotCMFSaved=true;//not saved for this time yet
		}
	}
	
	//Set end time for non fragmenting cases, no need to evolve for too long
	if(GammaLimitReached==false && gamma>=NONFRAGMENTGAMMALIMIT){
		GammaLimitReached=true;//setting flag
		Tmax=Min_2_double(t+1.0, Tmax); //evolution for one more dynamical time is allowed	
	}
	
	t+=dt; //increase time
	tcurrent+=dt*timeconversionfactor; //get real time

 }//end time evolution while loop

//Check if it totally collapsed or too much time elapsed or too small fragment to matter
if( (y0[0]<=ymin || t>Tmax || Mass<MinMass) && (Mass>IgnoreMass) && (HasFragmengted==false || IsSelfGravitating(Mass, TMean, R, MachedgeSq, gamma, LenientRhoFactor)) ){
	//Tell the Master that the cloud formed a star
	#if (ExtraVerboseOutput==true)
	cout<<rank<<":\t"<<"Slave finished with a star: Mass="<<Mass<<"\t Pos=["<<X0<<", "<<Y0<<", "<<Z0<<"]\t Vel=["<<Vx<<", "<<Vy<<", "<<Vz<<"] \t t="<<tcurrent<<" Myr \n";
	#endif
	message=FinishStarMessage; //message= "I made a star for you"
	mpierror=MPI_Ssend( &message, 1, MPI_INT, 0, SendMessageTag, MPI_COMM_WORLD);//send message
	//Get star parameters
	mpierror=MPI_Ssend( &Mass, 1, MPI_DOUBLE, 0, SendMassTag, MPI_COMM_WORLD); //mass
	mpierror=MPI_Ssend( &X0, 1, MPI_DOUBLE, 0, SendXTag, MPI_COMM_WORLD); //X
	mpierror=MPI_Ssend( &Y0, 1, MPI_DOUBLE, 0, SendYTag, MPI_COMM_WORLD); //Y
	mpierror=MPI_Ssend( &Z0, 1, MPI_DOUBLE, 0, SendZTag, MPI_COMM_WORLD); //Z
	mpierror=MPI_Ssend( &Vx, 1, MPI_DOUBLE, 0, SendVXTag, MPI_COMM_WORLD); //VX
	mpierror=MPI_Ssend( &Vy, 1, MPI_DOUBLE, 0, SendVYTag, MPI_COMM_WORLD); //VY
	mpierror=MPI_Ssend( &Vz, 1, MPI_DOUBLE, 0, SendVZTag, MPI_COMM_WORLD); //VZ
	mpierror=MPI_Ssend( &tcurrent, 1, MPI_DOUBLE, 0, SendTimeTag, MPI_COMM_WORLD); //time
	#if (VerboseOutput==true)
	cout<<rank<<":\t"<<"Slave sent a star\n";
	#endif
    } // end totally collapsed if 
     //It did not collapse
     else{
	//Tell the Master that the cloud did not form a star
	#if (VerboseOutput==true)
	cout<<rank<<":\t"<<"Slave finished without star\n";
	#endif
	message=FinishNoStarMessage; //message= "No star from this cloud"
	mpierror=MPI_Ssend( &message, 1, MPI_INT, 0, SendMessageTag, MPI_COMM_WORLD);//send message
	#if (VerboseOutput==true)
	cout<<rank<<":\t"<<"Slave sent message about not finding star\n";
	#endif
    }// end else, for case it did not collapse      
//Free variables 
free(outputstruct);
free(RelSizeVector);
free(Mach);
free(delta);
free(S);
free(rho_old);
free(T_old);
free(transform);
free(deltaK);
free(rho_T_field);
free(pointX);free(pointY);free(pointZ);
free(rho_field_new);
free(T_field_new);
free(y0);free(y1);free(y2);
free(Lvect);free(SqrtLvect);free(InvLvectSq);
free(NSpheres); free(SphereVolume);
free(renormstep); free(decaystep);
free(Xindex);free(Yindex);free(Zindex);
free(Xcoord);free(Ycoord);free(Zcoord);
free(MinusKIndex);
free(Klength);free(InvKRelLength);free(OneOverKdensity4PiKcube);free(SqrtKRelLength); 
#if (VerboseOutput==true)   
cout<<rank<<":\t"<<"Evolvefield cleanup finished ...\n";
#endif
        
}//end of code



//Calculates sound speed
inline double SoundSpeedSquare(double TMean, double gamma)
{
  return (gamma*kb/mu*TMean);            
       }
  
 //Critical density using Jeans criterion     
inline double RhoCrit(double TMean, double R, double MachedgeSq, double gamma){
	return (SoundSpeedSquare(TMean, gamma)*(1.0+MachedgeSq)/(R*R*FourPiOver3*Gconst*QVirial));
}

//Check if the cloud is self gravitating, using Jeans condition
inline bool IsSelfGravitating(double SphereMass, double TMean, double R, double MachedgeSq, double gamma, double StrictnessFactor){
    if( SphereMass >  StrictnessFactor*(FourPiOver3*pow(R,3.0)*RhoCrit(TMean,R,MachedgeSq, gamma)) ){ 
         return true;
         }
    else{
         return false;
         }
}

double FindCMFTime(double tcurrent)
{
	int j;
	double res;
	const double CMFTimes[CMFTIMENUM] = CMFTIMEVALUES; //times at which CMF and IMF are saved. 
	//Find nearest CMF observation time
	res=-1.0;//init
	if(tcurrent<=CMFTimes[0]){
		res=CMFTimes[0];
	}
	else{
		for ( j = 1; j < CMFTIMENUM; j++ ){
			if(CMFTimes[j-1]<tcurrent && CMFTimes[j]>=tcurrent){ //between these values
				res=CMFTimes[j];
			}
		}
	}

	return res;
}

long FindCMFTimeIndex(double tcurrent)
{
	int j;
	long res;
	const double CMFTimes[CMFTIMENUM] = CMFTIMEVALUES; //times at which CMF and IMF are saved. 
	//Find nearest CMF observation time
	res=-1;//init
	if(tcurrent<=CMFTimes[0]){
		res=0;
	}
	else{
		for ( j = 1; j < CMFTIMENUM; j++ ){
			if(CMFTimes[j-1]<tcurrent && CMFTimes[j]>=tcurrent){ //between these values
				res=j;
			}
		}
	}

	return res;
}

double SetOptimaldt(double dt,double renormstepMax,double dDeltaKMax){
	double newdt; // new value for dt/t0 to use
	double newrenormstepMax; //the value of renormstepMax after changing to newdt
	
	// First optimize for density variance steps (this is quadratic)
    newdt=dt*(dDeltaKOptimal/dDeltaKMax)*(dDeltaKOptimal/dDeltaKMax);
    newrenormstepMax=renormstepMax*newdt/dt; //linear relation
    // Check if this is acceptable for the dt/tau steps
	if (newrenormstepMax>RenormstepOptimal){
    	newdt=newdt*(RenormstepOptimal/newrenormstepMax); //linear dependence
	}
	
return newdt;	
}


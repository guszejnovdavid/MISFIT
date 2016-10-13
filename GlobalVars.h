/**************************************************************************************/
								/*  MISFIT Settings */
								
// This is where the parameters of MISFIT are set.

#ifndef GlobalVars_INCLUDED
#define GlobalVars_INCLUDED

/**************************************************************************/
/* Main Code Inputs */

//Maximum runtime in seconds (Zwicky has 48 hours walltime limit)
#define Max_RunTime 48*3500

//Sonic length in pc (sets the strength of turbulence: e.g. linewidth-size relation). Default is 0.1 pc (usual scale of turbulent filaments)
#define Rs 0.1
//Sonic mass in Msolar (sets T, provided Rs is given). Default T=10 K for Rs=0.1 pc means Ms=1.6
#define Ms 1.6
//Initial GMC mass (in Msolar)
#define StartMass 1.0E3

/**************************************************************************/
/* Physics Settings */

//Sets EOS of the gas to isothermal (only to compression, heating is handled differently), Default: true
#define IsothermalEOS true

//Enable protostellar heating based on Krumholz 2011, Default: true
#define HeatingAllowed true
//Set heating temperature to scale as R^{-0.5}, thus P=(...)T^4=(...)/R^2, default: true
#define R05Heating false
//Stochastic heating parameters (the heating gets a randomly (lognormal) assigned efficiency for each cloud)
#define StochasticHeating false //is it allowed?, default: false
#define HeatingSVal 1.0 //variance of the lognormal distribution, default: 1.0 (no physical reason)

/**************************************************************************/
/* Termination Conditions for Fragmentation Cascade */

//Stop size limit (evolution stops when a cloud reaches this relative size), as angular momentum support kicks in. 0.01 value based on Burkert_Bodenheimer_2000 average beta
#define MINIMUMSIZE 1.0E-2

// Minimum fragment mass allowed, fragments below are ignored, this is set by the opcaity limit from Low-Lynden-Bell 1976 to 0.007 Msolar
#define IgnoreMass 0.007
// Minimum fragment mass to consider for recursion, objects below this are assumed to collapse without fragmenting (default: same as IgnoreMass)
#define MinMass 0.007


/**************************************************************************/
/* Numerical Settings */

//Number of grid points per direction, powers of 2 are optimal for FFT,  USE POWER OF 2 FOR Ngrid!!!! This is the "resolution" of the simulation, runtime scales as Ngrid^3
#define Ngrid 32 //IMF usually converges at N=32

//number of size fragment scales to use (nmuber of sizes for MC spheres), default: 10
#define NScales 10

// Time used to initialize the turbulence (in frefall time), default: 2
#define InitTime 2.0
//Default time step used, in frefall time, default: 0.01
#define DEFAULTTIMESTEP 0.01
//Init time step (the first few time steps are very small so that the hierarchy for newly initialized scales can be established), in frefall time, default: 1E-6
#define FIRSTTIMESTEP 1E-6
//Number of small initializing steps these use stepsize FIRSTTIMESTEP, default: 3
#define FIRSTTIMESTEPNUM 3
// Optimal step sizes (<<1 due to linear approximation and better convergence to central limit theorem)
#define dDeltaKOptimal 0.10 //optimal step size for log density modes, default 0.10
#define RenormstepOptimal 0.10 //optimal step size for dt/tau, default 0.10


/**************************************************************************/
/* General Flags */


/*********************/
/* Debugging/verbose output  */

//Verbose output, maily for debugging
#define VerboseOutput false
// Very specific messages, mainly for debugging
#define ExtraVerboseOutput false
//Fix Random Seed, seed for all evolution calculations start at 0, FOR DEBUGGING, default=false
#define RandomSeedFixed false


/*********************/
/* Saving the evolution of self gravitating clouds  */

// Are we interested in the CMF evolution (if set to 1 the CMF is saved at predefined times) (Default true)
#define CMFEvolSave true
//CMF observation times (Myr)
#define CMFTIMEVALUES {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0}
//length of above vector (Must match the length of the vector above)
#define CMFTIMENUM 16
// Are we interested in the density fields? (Default false)
#define DensitySave true


/*********************/
/* Memory management  */

// Number of new lines allocated for StarMatrix and CoreMatrix at a given time (default: 1000)
#define LineIncreaseUnit 1000
// maximum number of doubles used by star and core matrices, before they get dumped to the hard drive
#define MaxAllocated 1E6
//Number of new structs added to InputStorage at a given time, default 50
#define StructIncreaseUnit 50
//Maximum number of cores to use (one computing clusters memory/core is fixed, with this one can limit memory usage).
//default: 0 (means no limit), the minimum is 2 (1 master, 1 slave) 
#define MaxCores 0


/*********************/
/* Extra Flags  */

//Is turbulence initialization allowed? (Default 1=YES)
#define InitAllowed 1

//Can clouds fragment (Default: true)
#define FragmentAllowed true

//Should the global parameters of the evolving parent be saved? Used to just see how the global parameters of a collapsing cloud evolve. (Default: false)
#define SaveParentEvol true

//When looking for self-gravitating substructures, should we start at large scales? True means that the search starts with the largest MC spheres, default: true
#define StartWithLargeFragments true

//Should skip the check that requires a fragment to be still self-gravitating when its heating is turned on, Default: true
#define SkipExtraCollapseCheck true

/**************************************************************************/
/* Constants (more or less motivated) */

//Turbulent mode strength b^2: S=log(1+b^2*Mach^2)
#define BSq 0.0625

//Virial parameter (kept constant)
#define QVirial 1.2 // standard value, using Jeans criterion if k~(Pi/2)/R. Using k~1/R yields Q=3 (Phil's paper uses that), but I think Q=1.2 is more realistic

//Limit for polytropic index, above which fragmentation is unlikely, so the evolution time is shortened. For a sphere 4/3 is the stability limiting case, that's the default.
#define NONFRAGMENTGAMMALIMIT 1.33

//Maximum relative mass of fragment to parent (does not make sense to have a "fragment" with 90% of the mass)
#define MaxFragmentRelativeMass 0.7

//Fraction of the critical mass required for collapse in STRICT mode (how high the density of a subregion should be before we allow it to fragment, relative to critical density)
#define StrictRhoFactor 1.01
//Fraction of the critical mass required for collapse in LENIENT (not STRICT) mode (e.g. checking at the end for a cloud if it is self-gravitating)
#define LenientRhoFactor 0.90 //this should ensure that if a core loses a small portion of its mass to a fragment which makes just below the thermal support limit, it still forms a star

const long NgridTotal=Ngrid*Ngrid*Ngrid; // Total number of grid points, used a LOT

/**************************************************************************/
/* Physical constants  */

//Grav. constant in m^2/s^2/Msolar*pc
#define Gconst 4325.69
//pc/Myr, used for conversiom
#define PcPerMyr 978.462
// Constant value for 4 Pi G in units of m^2/s^2/Msolar*pc
#define FourPiG 54358.2
//Time unit (t0=Sqrt(pc^3/Msolar/G)/Myr)
#define TimeUnit 14.877
//Boltzman constant in SI
#define kb 1.3806488E-23
//Mass of gas molecules (H2) in kg
#define mu 3.34524E-27 //2.0*(1.67262178E-27)

/**************************************************************************/
/* MPI communication  */

//MPI communication messages
	//From Slave
#define FinishStarMessage 100
#define FinishNoStarMessage 101
#define CloudFoundMessage 102
#define RecursionMessage 103
#define FieldsMessage 104
	//From Master
#define NewRunMessage 104
#define EndRunMessage 105
//MPI communication tags
#define SendInputTag 1
#define SendMessageTag 2
#define SendFieldsTag 3
	//sending object data
#define SendMassTag 10
#define SendXTag 11
#define SendYTag 12
#define SendZTag 13
#define SendVXTag 14
#define SendVYTag 15
#define SendVZTag 16
#define SendTimeTag 17
#define SendTimeIndexTag 18
#define SendMachTag 19
#define SendTemperatureTag 20
#define SendRadiusTag 21
#define SendRhoTag 22
#define SendMassInitTag 23


#endif

/**************************************************************************************/
								/*  MISFIT Master process */
								
// Manages the calculation (rank=0 process is the Master), collects the results and distributes tasks.

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <stddef.h>
#include <cstdio>

//Make it run with MPI for cluster
# include "mpi.h"

//Include own stuff
#include "Master.h"
#include "GlobalVars.h"
#include "InputDataManage.h"
#include "SumScalarInRegion.h"
#include "FindParameter_Newton.h"
#include "vectorfunctions.h"
#include "EvolveField.h"
#include "StarsAndCores.h"
#include "Restart.h"
#include "FileManage.h"

using namespace std;
    

/******************************************************************************/
// Master process that distributes tasks between Slaves
int Master(long Ntries){

/*******************************/
	/*Control variables*/
	
ofstream myfile; //Output file handle
time_t start, EndTime,currtime; //Time
bool OutOfTime=false; // set to true if the run should finish and save values to disk
long j; //loop variables
long Niter=0; //loop variable, counts how many GMCs have been done
long Nlines,NStruct; //used fro restarting, first is the number of lines in Restart_Inputs.txt, second is the number of input structs it corresponds to
long indexToCall; //index of the next run to be called
char filename[100]; //string used for filenames
bool IsRestart=false; //flag, set to 1 if this is a restart
bool *IsIdle; //vector, shows which slaves are idle
int NIdle; //number of idle slaves
int n_proc; //number of CPUs
int proc_ind; //index of slave to be called
int mpierror; //integer for mPI error messages
MPI_Status status; //used for determining sender of message
int sender; //ID of slave sending the message
int message; //message code from slave
long containerindex; //index of InputStorage into which the new struct is written

/*******************************/	
	/*Parameter variables*/
	
//Initial conditions
double T_Init; //initial temperature
double R_Init; //Initial radius vector
double Mach_Init; //Initial edge Mach number
double rho_Init; //initial density
double* dummyvector; //dummy vector used for getting proper normalization
double normfactor;
struct inputstruct_s *inputstruct, *default_inputstruct;
#if (DensitySave==true)
	//structure used to send field data to master
	struct densitystruct_s *fieldstruct;
#endif
//Star and core data, used as output
double ObjectMass;
double ObjectX, ObjectY, ObjectZ;
double ObjectVX, ObjectVY, ObjectVZ;
double ObjectTime;
//extra for clouds
long ObjectTimeindex;
double ObjectMassInit,ObjectMach, ObjectTemperature, ObjectRadius, ObjectRhoAverage;



/*******************************/
	/*Shorthands*/
	
const double FourPiOver3=4.0/3.0*M_PI;

/***********************************************************************************************************************************************/
			/* Set Up */

// set up initial condition for one GMC

//if restart run
	// read restart control
	// read stored inputs for new runs
// if first run
	// initialize data storage with one GMC (using initial conditions)
	
/*******************************/	
/* Control variables */
#if (VerboseOutput==true)
cout<<"Master starting\n";
#endif

//get number of CPU	
MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
if(n_proc>MaxCores && MaxCores>0){ //check if we are limiting them
	n_proc=MaxCores;
}
//Init idle vector
IsIdle=vector_bool((long)n_proc);
IsIdle[0]=false; //Master is always working
for(j=1;j<(long)n_proc;j++){
	IsIdle[j]=true; //all slaves start out without a job
}
	
/*******************************/	
/* Set up initial conditions */

//Allocate memory
dummyvector=(double *) calloc(NgridTotal, sizeof (double)); //dummy vector of N^3 size
//input structures can be several MB large
inputstruct=(struct inputstruct_s*) calloc(1,sizeof(struct inputstruct_s));
 #if (DensitySave==true)
	 fieldstruct=(struct densitystruct_s*) calloc(1,sizeof(struct densitystruct_s)); //allocate memory for field structure used to send data to master
#endif
default_inputstruct=(struct inputstruct_s*) calloc(1,sizeof(struct inputstruct_s));

//Time set up
time(&start); //time of start
EndTime=start+Max_RunTime; //after this no new slave process will be started

//Assign values
T_Init=QVirial/2.0*Ms*Gconst*mu/(kb*Rs); //Initial temperature, set by the sonic mass and scale (~10 K), using gamma=1 (true for large cloud)
vector_set( dummyvector, 1.0, NgridTotal );//setting initial density as constant 1
normfactor=FourPiOver3/SumScalar_Sphere(dummyvector, Ngrid,(-1.0), 1.0, (-1.0), 1.0, (-1.0), 1.0, 1.0, 0.0, 0.0, 0.0); // get summation normalization factor
R_Init=FindR(StartMass); //Initial radius in pc
Mach_Init=sqrt(R_Init/Rs); // Mach scales as sqrt(R)
rho_Init=normfactor* StartMass/(FourPiOver3*pow(R_Init,3)); //Initial density

//Fill up initial struct
vector_set(default_inputstruct[0].rho_field, rho_Init, NgridTotal );
vector_set(default_inputstruct[0].T_field, T_Init, NgridTotal );
default_inputstruct[0].R0=R_Init;
default_inputstruct[0].X0=0.0;default_inputstruct[0].Y0=0.0;default_inputstruct[0].Z0=0.0; //center of reference frame
default_inputstruct[0].Machedge=Mach_Init;
default_inputstruct[0].RelativeVelocityX=0.0;default_inputstruct[0].RelativeVelocityY=0.0;default_inputstruct[0].RelativeVelocityZ=0.0; //at rest in reference frame
default_inputstruct[0].InitScale=InitAllowed*10.0; //initialize turbulence on all scales, if allowed
default_inputstruct[0].tparent=0.0; //starts from t=0
//Create output folder if it does not exist
CreateDir("Output");

//Check if the run has already finished
ContructFilePath("Finish.txt", "Output", filename);
if(file_exists (filename)){
	cout<<"Finish.txt found in Output folder. Run already finished??? Terminating.\n";
	//Tell all slaves to terminate
	for(proc_ind=1;proc_ind<n_proc;proc_ind++){// all slave processes
		message=EndRunMessage; //message= no more runs, terminate
		mpierror=MPI_Ssend( &message, 1, MPI_INT, proc_ind, SendMessageTag, MPI_COMM_WORLD);//send message
	}
	return 1; //terminate the run
}

//Check if this is a restart
ContructFilePath("RestartControl.txt", "Output", filename);
if(file_exists (filename)){
	IsRestart=true;
	cout<<"Restart control file found.\n";
}

if(IsRestart){
	
	/********************************************/
		/* Restart Run  */

	ReadRestartControl(&Niter); //Get Niter from contrl file
	ContructFilePath("Restart_Inputs.txt", "Output", filename); // path for stored input data
	Nlines=CountNumberOfLines(filename); //number of lines in file
	NStruct=Nlines/NDataPerStruct; //number of structures stored
	InitInputDataStorage(NStruct+1);// sets up InputStorage, adds an extra space for new input
	ReadRestartFile(NStruct); //reads input structs and stores them in InputStorage
}
else{

	/********************************************/
		/* First Run  */

	Niter=0; //first run
	InitInputDataStorage(2);// sets up InputStorage, adds an extra space for new input
	containerindex=FirstUnused(); //first empty storage, will be 0 in this case
	PutIntoInputStorage(default_inputstruct,containerindex); //put new GMC into storage, use first slot
}
InitStarsAndCores(); // Sets up star and core storage matrices

/***********************************************************************************************************************************************/
				/*  MASTER Loop  */
				
//Ntries loop start

	//start highest priority initial dataset
				
	//Repeat
		// Wait for message from any slave
		
		// if  "finished"
			//store star parameters
			// store slave ID in "idle" vector, increase NAvailableCores by 1
		// if "recursion"
			// receive input data struct from slave, increase NInputData
			//store input data struct, calculate priority and store it in vector
		// if "core found"
			//store core data
		
		//if NAvailableCores>=1 and NInputData>=1 and time<Endtime
			// call a slave with highest priority initial dataset
			
		// if NAvailableCores=NCores-1 (everyone idle, master always busy) 
			// terminate loop
			
	//EndRepeat
	
	// if time>EndTime
		// save to initial condition data torestart file
		// create restart control file
		// terminate Ntries loop
	//else
		// set up one GMC with initial parameters (we start the next GMC)
		
//End of Ntries loop
/*************************************************************************************/
#if (VerboseOutput==true)
cout<<"Master loop starting.\n";
#endif
//Loop through all GMCs
while(Niter<Ntries){
	
	//Find which parameters to use for first run
	indexToCall=NextToCall();
	//Find the a slave that can be called
	proc_ind=FindIdle(IsIdle, n_proc);
	#if (ExtraVerboseOutput==true)
	cout<<"First run call init.\n";
	#endif
	//Call first slave run, tell it to prepare for new data
	message=NewRunMessage; //message= prepare for new run
	mpierror=MPI_Ssend( &message, 1, MPI_INT, proc_ind, SendMessageTag, MPI_COMM_WORLD);//send message
	//send it the input data from storage
	#if (ExtraVerboseOutput==true)
	cout<<"Master sends first data.\n";
	#endif
	mpierror=MPI_Ssend( &InputStorage[indexToCall], 1, mpi_inputstruct_type, proc_ind, SendInputTag, MPI_COMM_WORLD);
	IsIdle[proc_ind]=false; //no longer idle
	IsUsed[indexToCall]=false; //we already used this
	NInputCurrent--;
	
	//Start Master control loop, only stops when broken
	while(1){
		//Wait for message back from any slave
		#if (VerboseOutput==true)
		cout<<"Master waits for message\n";
		#endif
		mpierror=MPI_Recv( &message, 1, MPI_INT, MPI_ANY_SOURCE, SendMessageTag, MPI_COMM_WORLD, &status);
		sender = status.MPI_SOURCE; //determine source
		
		//Slave finished, there is a star
		if(message==FinishStarMessage){
			#if (VerboseOutput==true)
			cout<<"Master is receiving a star from Slave "<< sender<<"\n";
			#endif
			//Check if there is enough space for one new dataline, if not, then allocate it
			if(NStarsCurrent==NStarsAllocated){
				IncreaseStarMatrix();
			}
			//Get star parameters
			mpierror=MPI_Recv( &ObjectMass, 1, MPI_DOUBLE, sender, SendMassTag, MPI_COMM_WORLD, &status); //mass
			mpierror=MPI_Recv( &ObjectX, 1, MPI_DOUBLE, sender, SendXTag, MPI_COMM_WORLD, &status); //X
			mpierror=MPI_Recv( &ObjectY, 1, MPI_DOUBLE, sender, SendYTag, MPI_COMM_WORLD, &status); //Y
			mpierror=MPI_Recv( &ObjectZ, 1, MPI_DOUBLE, sender, SendZTag, MPI_COMM_WORLD, &status); //Z
			mpierror=MPI_Recv( &ObjectVX, 1, MPI_DOUBLE, sender, SendVXTag, MPI_COMM_WORLD, &status); //VX
			mpierror=MPI_Recv( &ObjectVY, 1, MPI_DOUBLE, sender, SendVYTag, MPI_COMM_WORLD, &status); //VY
			mpierror=MPI_Recv( &ObjectVZ, 1, MPI_DOUBLE, sender, SendVZTag, MPI_COMM_WORLD, &status); //VZ
			mpierror=MPI_Recv( &ObjectTime, 1, MPI_DOUBLE, sender, SendTimeTag, MPI_COMM_WORLD, &status); //time
			if(ObjectMass>MinMass){
			//Store star parameters
			AssignStarMatrixValues(ObjectMass, ObjectX, ObjectY, ObjectZ, ObjectVX, ObjectVY, ObjectVZ, ObjectTime);
			}
			//Slave is now idle
			IsIdle[sender]=true;
			#if (ExtraVerboseOutput==true)
			cout<<"0:\t Master received a star with Mass="<<ObjectMass<<"\t Pos=["<<ObjectX<<", "<<ObjectY<<", "<<ObjectZ<<"]\t Vel=["<<ObjectVX<<", "<<ObjectVY<<", "<<ObjectVZ<<"] \t t="<<ObjectTime<<" Myr \n";	
			#endif
		}
		
		//Slave finished, there is no star
		if(message==FinishNoStarMessage){
			#if (VerboseOutput==true)
			cout<<"Master receives no star termination from Slave "<< sender<<"\n";
			#endif
			//Slave is now idle
			IsIdle[sender]=true;	
		}
		
		//Slave found a cloud that exists at one of the pre-defined observation times
		if(message==CloudFoundMessage){
			#if (VerboseOutput==true)
			cout<<"Master receives a cloud from Slave "<< sender<<"\n";
			#endif
			//Get time index
			mpierror=MPI_Recv( &ObjectTimeindex, 1, MPI_LONG, sender, SendTimeIndexTag, MPI_COMM_WORLD, &status);
		    //Check if there is enough space for one new dataline, if not, then allocate it
			if(NCoresCurrent[ObjectTimeindex]==NCoresAllocated[ObjectTimeindex]){
				IncreaseCoreMatrix(ObjectTimeindex);
			}
			//Get cloud parameters
			mpierror=MPI_Recv( &ObjectMass, 1, MPI_DOUBLE, sender, SendMassTag, MPI_COMM_WORLD, &status); //mass
			mpierror=MPI_Recv( &ObjectX, 1, MPI_DOUBLE, sender, SendXTag, MPI_COMM_WORLD, &status); //X
			mpierror=MPI_Recv( &ObjectY, 1, MPI_DOUBLE, sender, SendYTag, MPI_COMM_WORLD, &status); //Y
			mpierror=MPI_Recv( &ObjectZ, 1, MPI_DOUBLE, sender, SendZTag, MPI_COMM_WORLD, &status); //Z
			mpierror=MPI_Recv( &ObjectVX, 1, MPI_DOUBLE, sender, SendVXTag, MPI_COMM_WORLD, &status); //VX
			mpierror=MPI_Recv( &ObjectVY, 1, MPI_DOUBLE, sender, SendVYTag, MPI_COMM_WORLD, &status); //VY
			mpierror=MPI_Recv( &ObjectVZ, 1, MPI_DOUBLE, sender, SendVZTag, MPI_COMM_WORLD, &status); //VZ
			//extra parameters
			mpierror=MPI_Recv( &ObjectMassInit, 1, MPI_DOUBLE, sender, SendMassInitTag, MPI_COMM_WORLD, &status); //initial mass
			mpierror=MPI_Recv( &ObjectMach, 1, MPI_DOUBLE, sender, SendMachTag, MPI_COMM_WORLD, &status); //mach number
			mpierror=MPI_Recv( &ObjectTemperature, 1, MPI_DOUBLE, sender, SendTemperatureTag, MPI_COMM_WORLD, &status); //Temperature
			mpierror=MPI_Recv( &ObjectRhoAverage, 1, MPI_DOUBLE, sender, SendRhoTag, MPI_COMM_WORLD, &status); //average density
			mpierror=MPI_Recv( &ObjectRadius, 1, MPI_DOUBLE, sender, SendRadiusTag, MPI_COMM_WORLD, &status); //radius	

			//Store cloud parameters
			AssignCoreMatrixValues(ObjectTimeindex, ObjectMassInit, ObjectMass, ObjectX, ObjectY, ObjectZ, ObjectVX, ObjectVY, ObjectVZ, ObjectMach, ObjectTemperature, ObjectRadius, ObjectRhoAverage);	
		}
		
		//Slave found a cloud that exists at one of the pre-defined observation times
		if(message==FieldsMessage){
			#if (VerboseOutput==true)
			cout<<"Master receives field data from Slave "<< sender<<"\n";
			#endif
			//Get time index
			mpierror=MPI_Recv( &ObjectTimeindex, 1, MPI_LONG, sender, SendTimeIndexTag, MPI_COMM_WORLD, &status);
			#if (ExtraVerboseOutput==true)
			cout<<"Timeindex for fields "<< ObjectTimeindex<<" from Slave "<<sender<<"\n";
			#endif
		   	//Get data
			mpierror=MPI_Recv( fieldstruct, 1, mpi_densitystruct_type, sender, SendFieldsTag, MPI_COMM_WORLD, &status);
			#if (ExtraVerboseOutput==true)
			cout<<"Master successfuly received fields data from Slave "<< sender<<"\n";
			#endif
			// Save data to file
			SaveFields(fieldstruct,ObjectTimeindex);
			#if (VerboseOutput==true)
			cout<<"Fields saved \n";
			#endif
		}
		
		//Slave found self-gravitatinng substructure, sends parameter for new run
		if(message==RecursionMessage){
			#if (VerboseOutput==true)
			cout<<"Master receives recursion data from Slave "<< sender<<"\n";
			#endif
			//Check if we have enough space to take one more
			if(NInputAllocated==NInputCurrent){
				IncreaseInputStorage(); //allocate more memory
			}
			//Find first empty container
			containerindex=FirstUnused();
			//Get data
			mpierror=MPI_Recv( inputstruct, 1, mpi_inputstruct_type, sender, SendInputTag, MPI_COMM_WORLD, &status);
			PutIntoInputStorage(inputstruct,containerindex); //put new cloud data into storage
			#if (ExtraVerboseOutput==true)
			cout<<"Master stored recursion data in InputStorage\n";
			#endif
		}
			
		/**********************************************/	
		/*Check up on core and star matrices*/
		if(message==CloudFoundMessage || message==FinishStarMessage){
			//Check if we are over our limit for the core and star matrices
		    if(TooLargeStarsAndCoresMatrices()){
		    	#if (VerboseOutput==true)
		    	cout<<"Too many star/clouds. Saving data to files to save memory.\n";
		    	#endif
		    	SaveAndFreeStarsAndCoresMatrices(false); //save to file and free memory
		    	InitStarsAndCores(); //reinitialize matrices
			}
		}
		
		//Get the number of idle slaves
		NIdle=SumIdle(IsIdle, n_proc);
		
		/************************************************/
		/* Call new slaves, if there are idle ones  */
		//Find which parameters to use for first run
		indexToCall=NextToCall();
		#if (ExtraVerboseOutput==true)
		cout<<"Index of data struct to call: "<<indexToCall<<"\t Number of idle slaves: "<<NIdle<<"\n";
		#endif
		if((indexToCall!=-1)&&(NIdle>0)&&(OutOfTime==false)){ //there is an idle slave, there is something to be calculated and we have time
			//Find the a slave that can be called
			proc_ind=FindIdle(IsIdle, n_proc);
			//Call a slave run, tell it to prepare for new data
			message=NewRunMessage; //message= Wake up! prepare for new run
			mpierror=MPI_Ssend( &message, 1, MPI_INT, proc_ind, SendMessageTag, MPI_COMM_WORLD);//send message
			//Put the slave to work!
			mpierror=MPI_Ssend( &InputStorage[indexToCall], 1, mpi_inputstruct_type, proc_ind, SendInputTag, MPI_COMM_WORLD);
			IsIdle[proc_ind]=false; //no longer idle
			IsUsed[indexToCall]=false; //we already used this
			NInputCurrent--;
			#if (VerboseOutput==true)
			cout<<"Master sent recursion data to slave "<<proc_ind<<"\n";
			#endif
		}
		
		/***********************************************/
		// Look at the time!
		time(&currtime); //time now
		if(currtime>EndTime){
			cout<<"Master ran out of time, no more new runs started.\n";
			OutOfTime=true; // start to finish up
		}
		
		
		/**********************************************/
		// Terminate if every slave is done and there are nothing left to run
		NIdle=SumIdle(IsIdle, n_proc);
		#if (VerboseOutput==true)
		cout<<"Number of idle slaves: "<<NIdle<<"\t Used input storage: "<<NInputCurrent<<"/"<<NInputAllocated<<"\n";
		#endif
		if((NIdle==n_proc-1)&&(NInputCurrent==0)){
			#if (VerboseOutput==true)
			cout<<"Master: all recursions done, breaking loop\n";
			#endif
			break; //break loop
		}
		
		/**********************************************/
		// Terminate if out of time
		if(OutOfTime && (NIdle==n_proc-1)){ //ran out of time, but everyone finished
			cout<<"Ran out of time, terminating...\n";
			break;	
		}
	}
	
	//if we ran out of time, then we need to save our data and end this
	if(OutOfTime){
		#if (VerboseOutput==true)
		cout<<"Master ran out of time, saving data started...\n";
		#endif
		SaveAndFreeStarsAndCoresMatrices(false); //save to file and free memory, not the end of data for this GMC
		#if (VerboseOutput==true)
		cout<<"Stars and clouds saved.\n";
		#endif
		//save input structures
		SaveAndFreeInputStorage(); //save and set free...
		#if (VerboseOutput==true)
		cout<<"Input storage saved\n";
		#endif
		//save current position in loop
		SaveRestartControl(Niter);
		#if (VerboseOutput==true)
		cout<<"Restart control saved, breaking loop\n";
		#endif
		//break loop
		break;		
	}
	else{
		cout<<"Master finished for GMC "<< Niter<<"\n";
		//Save results to file
		SaveAndFreeStarsAndCoresMatrices(true); //save to file and free memory, is the end of data for this GMC
		InitStarsAndCores(); //reinitialize matrices
		
		Niter++;//go for next GMC
		//Put a new initial GMC into the input storage
		containerindex=FirstUnused(); //first empty storage, will be 0 in this case
		PutIntoInputStorage(default_inputstruct,containerindex); //put new GMC into storage, use first slot
	}	
}
#if (VerboseOutput==true)	
cout<<"Master loop ended\n";
#endif

/***********************************************************************************************************************************************/
				/*  Finishing up  */

//Tell all slaves to terminate
for(proc_ind=1;proc_ind<n_proc;proc_ind++){// all slave processes
	message=EndRunMessage; //message= no more runs, terminate
	mpierror=MPI_Ssend( &message, 1, MPI_INT, proc_ind, SendMessageTag, MPI_COMM_WORLD);//send message
}
#if (VerboseOutput==true)
cout<<"Slaves put to sleep\n";
#endif

if(OutOfTime){
	cout<<"Ran out of time. Data saved.\n";
	}
else{
		//Create end of run file
		ContructFilePath("Finish.txt", "Output", filename);
		WriteLineToTextFile("DONE", filename);
		//Remove inout files if they exist
		ContructFilePath("RestartControl.txt", "Output", filename);
		if(file_exists (filename)){
			cout<<"Removing old restart control file...\n";
			remove(filename);
		}

		
	}
cout<<"MISFIT run finished, clean up started...\n";

// Clean up
free(IsIdle);
//free(filename);
free(dummyvector);
free(inputstruct);
free(default_inputstruct);

//time at end
time(&currtime); //time now
cout<<"MISFIT calculation took "<< difftime(currtime, start) <<" seconds. \n";
cout<<"Exiting...\n";
return 0;
//End of code
}


//Finds the first idle CPU
int FindIdle(bool* IsIdle, int n_proc){
	int i;
	
	for(i=1;i<n_proc;i++){// all slave processes
		if(IsIdle[i]){
			return i; //index of first idle
		}
	}
	//cout<<"No idle slave.";
	return 0; // no idle process was found, return zero
}

//Finds the fnumber of idle slaves
int SumIdle(bool* IsIdle, int n_proc){
	int i,tot=0;
	
	for(i=1;i<n_proc;i++){// all slave processes
		if(IsIdle[i]){
			tot++;
		}
	}
	return tot; 
}
 
 
 

/**************************************************************************************/
			/*  Main Routine for Minimalist Star Formation with Intrinsic Turbulence (MISFIT) */
			
// This starts the calculation, initialized MPI and calls the Master and Slave processes.

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

//Make it run with MPI for cluster
# include "mpi.h"

//Include own stuff
#include "GlobalVars.h"
#include "InputDataManage.h"
#include "Master.h"
#include "Slave.h"
       

/******************************************************************************/
// Caller routine, makes it easy to switch to MPI version
int main(int argc, char *argv[]) {

/******************************************************************************/


using namespace std;
ofstream myfile; //Output file handle

/* Constants and parameters*/
long Ntries; //number of tries for each mass
long index; //index of the process, useful for multi_runs

int returnval; //returned value


/******************************************************************/
     /* Check input */
if ( argc != 2 ) // 3 arguments is the default (Ntries, Nmass, index)
{
    return 0;
}
else 
{
Ntries = strtol(argv[1],NULL,0); //setting endpoint of the calculation
}

cout<<"MISFIT Starting...\n";

#if (VerboseOutput==true)
cout<<"MISFIT MPI Init.\n";
#endif

//  Initialize MPI.
MPI::Init ( argc, argv );

//  Get the individual process ID.
index = MPI::COMM_WORLD.Get_rank();

#if (VerboseOutput==true)
cout<<"MPI datatypes init started.\n";
#endif

//Make Input struct passable by MPI
InitInputDataType();
//Make density struct passable by MPI
InitDensityStructDataType();

#if (VerboseOutput==true)
cout<<index<<": MPI datatypes initialized.\n";
#endif

if(index==0){
	//index=0 is for MASTER process that will distribute tasks between SLAVEs
	cout<<"Master called with rank:"<<index<<"\n";
	returnval=Master(Ntries);
}
else
{
	if(index<MaxCores || MaxCores==0){
	//index!=0 for SLAVE process that will do the actual calculation
		cout<<"Slave called with rank:"<<index<<"\n";
		returnval=Slave();
	}
	else{
		cout<<index<<": Terminate to conserve memory \n";
		returnval=0;
	}

}

//Free MPI datatype
MPI_Type_free(&mpi_inputstruct_type);

//free vectors
//free(offsets);free(blocklengths);

#if (VerboseOutput==true)
cout<<index<<": MPI Finalize.\n";
#endif

//Finalize
MPI::Finalize ( );
	

return returnval;

}

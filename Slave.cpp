/**************************************************************************************/
								/*  MISFIT Slave process */
								
// This process (rank!=0) uses the inputs received from Master to calculate the evolution of a single cloud.

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
#include "InputDataManage.h"
#include "EvolveField.h"


using namespace std;

//Include global variables
#include "GlobalVars.h"
    

/******************************************************************************/
// Slave routine, it calls a cloud calculation run and waits for new instructions
int Slave(){
	
	using namespace std;
	
	
	struct inputstruct_s *inputstruct; // to store input data
	int i; //loop
	int message; //message code 
	MPI_Status status; //MPI message status
	int mpierror; //integer for mPI error messages
	int rank; // process rank
	
// Slave Loop starts
	
	//Wait for message from Master
	
	//if new run
		//call EvolveField
	
	//if end run
		// break loop
		
//End of loop
	

/******************************************************************************/
//  Get the individual process ID.
rank = MPI::COMM_WORLD.Get_rank();
#if (VerboseOutput==true)
cout<<rank<<":\t"<<"Slave loop starting.\n";
#endif
//Allocate memory
inputstruct=(struct inputstruct_s*) calloc(1,sizeof(struct inputstruct_s));

//Initialize random number generator
#if (RandomSeedFixed==false)
srand(time(NULL)*(1+rank));
for(i=0;i<10000;i++){
	rand(); //throw away the first few
}
#endif


	while(1)
	{
		#if (VerboseOutput==true)
		cout<<rank<<":\t"<<"Slave waiting for message\n";
		#endif
		//Wait for message back from Master (index 0)
		mpierror=MPI_Recv( &message, 1, MPI_INT, 0, SendMessageTag, MPI_COMM_WORLD, &status);
		#if (VerboseOutput==true)
		cout<<rank<<":\t"<<"Slave received message\n";
		#endif
		if(message==NewRunMessage){
		mpierror=MPI_Recv(inputstruct, 1, mpi_inputstruct_type, 0, SendInputTag, MPI_COMM_WORLD, &status);
		#if (VerboseOutput==true)
		cout<<rank<<":\t"<<"Slave received data\n";
		#endif
		//Calling time evolving routine
		EvolveField(inputstruct);
		}
		
		if(message==EndRunMessage){
		#if (VerboseOutput==true)
		cout<<rank<<":\t"<<"Slave breaking loop\n";
		#endif
		//Break loop
		break;
		}	
	}
free(inputstruct);
	
return 0;
}

#ifndef InputDataManage_INCLUDED
#define InputDataManage_INCLUDED

#include "GlobalVars.h"
#include <stdio.h> //for having  file variables
# include "mpi.h"

// Declaring the input parameter struct type and creating the MPI passable version
typedef struct inputstruct_s {
        double rho_field[NgridTotal];
        double T_field[NgridTotal];
        double R0;
        double X0;double Y0;double Z0;
        double Machedge;
        double RelativeVelocityX; double RelativeVelocityY; double RelativeVelocityZ;
		double InitScale;
		double tparent;
		};  
		
// Declaring the density field struct type and creating the MPI passable version
typedef struct densitystruct_s {
        double rho_field[NgridTotal];
        double T_field[NgridTotal];
        double xcoord[NgridTotal];
        double ycoord[NgridTotal];
        double zcoord[NgridTotal];
        double volume;
		};  


#ifndef InputDataManage_DEFINED //Only use them if we are not definiing them locally
	/*****************************************/
			//Global variables
			
	extern MPI_Datatype mpi_inputstruct_type; //declare the new mpi type of the input structure 
	extern MPI_Datatype mpi_densitystruct_type; //declare the new mpi type of the density structure 
			
	// Input storage data
	extern struct  inputstruct_s  *InputStorage; //vector, containing the input parameter stuctures used for calling new runs.
	extern double *InputPriority; //vector containing the priority of each run
	extern bool *IsUsed; // vector showing which allocated spaces are used
	extern long NInputAllocated; //Number of structures allocated in InputStorage
	extern long NInputCurrent; //Number of structures stored in InputStorage
	
	extern long NDataPerStruct;
	

#endif

using namespace std;

/**********************************************/
	//Functions to manage Input Data Storage
	

		//Changing size
			
//Increases the number of structs by StructIncreaseUnit
void IncreaseInputStorage();

//Free Input Storage, used when we are done calculating
void FreeInputStorage();

//Save Input data and free memory, used when saving for later restart
void SaveAndFreeInputStorage();

		//Saving to disk
				
//Save input parameters for later restart
void SaveInputStorage(char* filename);


		//Storing data
		
//Store input stucture	
void PutIntoInputStorage(struct inputstruct_s* inputstruct, long pos_index);


//Make input data struct passabel by MPI	
void InitInputDataType();

//Make density data struct passabel by MPI	
void InitDensityStructDataType();

//Initialize data stoarge
void InitInputDataStorage(long Nitem);

//Calculate priority for input dataset
double PriorityFunction(struct inputstruct_s* inputstruct);

//Finds the index of the dataset with the highest priority
long NextToCall();

//Finds first unused slot in the dataset
long FirstUnused();

#endif

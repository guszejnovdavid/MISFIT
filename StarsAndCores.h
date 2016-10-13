#ifndef StarsAndCores_INCLUDED
#define StarsAndCores_INCLUDED

#include <stdio.h> //for having  file variables
#include "InputDataManage.h"

#ifndef StarsAndCores_DEFINED //Only use them if we are not definiing them locally
	/*****************************************/
			//Global variables
			
	// Data for stars
	extern double** StarMatrix; //matrix containing the parameters of formed protostars (size: 8 x Nstars).
						   // Mass  X0  Y0  Z0  RelativeVelocityX  RelativeVelocityY  RelativeVelocityZ  tcurrent
	extern long NStarsAllocated; //Number of lines allocated for StarMatrix
	extern long NStarsCurrent; //Number of lines filled out in StarMatrix
	
	// Core data
	extern double*** CoreMatrix; //matrix containing the parameters of formed protostars (size: CMFTIMENUM x 8 x Ncores).
						   // Mass  X0  Y0  Z0  RelativeVelocityX  RelativeVelocityY  RelativeVelocityZ  tcurrent
	extern long *NCoresAllocated; //Number of lines allocated in the core matrix for different times
	extern long *NCoresCurrent; //Number of lines filled out in the core matrix for different times
	
	//Restart file pointers
	extern FILE *restartfile; // pointer to the restart file
#endif

using namespace std;

/**********************************************/
	//Functions to manage global variables
	

		//Changing size
			
//Increases the number of lines in StarMatrix by LineIncreaseUnit
void IncreaseStarMatrix();

//Free StarMatrix
void FreeStarMatrix();

//Increases the number of lines in CoreMatrix by LineIncreaseUnit for a specific time
void IncreaseCoreMatrix(long timeindex);

//Free CoreMatrix
void FreeCoreMatrix();

//Save matrix data and free memory
void SaveAndFreeStarsAndCoresMatrices(bool IsEnd);

//Check if memory limit for star and core matrices is exceeded
bool TooLargeStarsAndCoresMatrices();

		//Saving to disk
				
//Save StarMatrix
void SaveStarMatrix(char* filename, bool IsEnd);

//Save StarMatrix
void SaveCoreMatrix(long timeindex, char* filename, bool IsEnd);

//Save field data
void SaveFields(struct densitystruct_s* fieldstruct, long ObjectTimeindex);

		//Storing data
		
//Store parameters of a protostar		
void AssignStarMatrixValues(double Mass, double X0, double Y0, double Z0, double RelativeVelocityX, double RelativeVelocityY, double RelativeVelocityZ, double tcurrent);

//Store parameters of a core		
void AssignCoreMatrixValues(long timeindex, double MassInit, double Mass, double X0, double Y0, double Z0, double RelativeVelocityX, double RelativeVelocityY, double RelativeVelocityZ, double Machedge, double Temperature, double Radius, double RhoAverage);

//Init	
void InitStarsAndCores();

//Naming of core/cloud data files
void CMFFileName(char* fullname, char* filename, double CMFTime);
#endif

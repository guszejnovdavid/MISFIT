/**************************************************************************************/
			/*  Handling star and cloud output data */
			
// Collection of routines used by Master for storing output data

#include <stdlib.h>
#include "vectorfunctions.h"
#define StarsAndCores_DEFINED // tells StarsAndCores.h not to define the external variables as we will do it here
#include "GlobalVars.h"
#include "StarsAndCores.h"
#include "FileManage.h"
#include <fstream>
#include <iostream>

/*****************************************/
		//Global variables
		
// Data for stars
double** StarMatrix; //matrix containing the parameters of formed protostars (size: 8 x Nstars).
					   // Mass  X0  Y0  Z0  RelativeVelocityX  RelativeVelocityY  RelativeVelocityZ  tcurrent
long NStarsAllocated; //Number of lines allocated for StarMatrix
long NStarsCurrent; //Number of lines filled out in StarMatrix

// Core data
double*** CoreMatrix; //matrix containing the parameters of formed protostars (size: CMFTIMENUM x 8 x Ncores).
					   // Mass  X0  Y0  Z0  RelativeVelocityX  RelativeVelocityY  RelativeVelocityZ  tcurrent
long *NCoresAllocated; //Number of lines allocated in the core matrix for different times
long *NCoresCurrent; //Number of lines filled out in the core matrix for different times

//Restart file  pointers
FILE *restartfile=NULL; // pointer to the restart file

const double CMFTimes[CMFTIMENUM] = CMFTIMEVALUES; //times at which CMF and IMF are saved.


/*************************************************************/
		//Changing size
		
		
//Increases the number of lines in StarMatrix by LineIncreaseUnit
void IncreaseStarMatrix(){
	long i;
	double *temppointer;
	
	for(i=0;i<8;i++){ //reallocating memory for all lines, increasing the number of elements
		temppointer=(double*) realloc(StarMatrix[i], (NStarsAllocated+LineIncreaseUnit)* sizeof(double));
		if(temppointer!=NULL){ //successful allocation
			StarMatrix[i]=temppointer;
		}
		else
		{
			cout<<"Memory allocation failure! Exiting!";
			exit(EXIT_FAILURE);
		}
	}
	NStarsAllocated+=LineIncreaseUnit; //note the number of elements	
}
//Free StarMatrix
void FreeStarMatrix(){
	long i;
	
	for(i=0;i<8;i++){ 
		free(StarMatrix[i]);
	}
	free(StarMatrix);
	StarMatrix=NULL; //none assigned
	NStarsAllocated=0; //zero size
	NStarsCurrent=0; //zero stars
}

//Increases the number of lines in CoreMatrix by LineIncreaseUnit for a specific time
void IncreaseCoreMatrix(long timeindex){
	long i;
	double *temppointer;
	
	for(i=0;i<12;i++){ //reallocating memory for all lines, increasing the number of elements
		temppointer=(double*) realloc(CoreMatrix[timeindex][i], (NCoresAllocated[timeindex]+LineIncreaseUnit)* sizeof(double));
		if(temppointer!=NULL){ //successful allocation
			CoreMatrix[timeindex][i]=temppointer;
		}
		else
		{
			cout<<"Memory allocation failure! Exiting!";
			exit(EXIT_FAILURE);
		}
	}
	NCoresAllocated[timeindex]+=LineIncreaseUnit; //note the number of elements	
}
//Free CoreMatrix
void FreeCoreMatrix(){
	long i,j;
	
	for(i=0;i<CMFTIMENUM;i++){ 
		for(j=0;j<12;j++){ 
			free(CoreMatrix[i][j]);
		}
		free(CoreMatrix[i]);
		NCoresAllocated[i]=0;
		NCoresCurrent[i]=0;	
	}
	free(NCoresAllocated);
	free(NCoresCurrent);
	free(CoreMatrix);
	CoreMatrix=NULL; //none assigned
}

//Save matrix data and free memory
void SaveAndFreeStarsAndCoresMatrices(bool IsEnd){
	long l;
	char filename[100];
	char Corefilename[100]; //file containing the cores at a chosen time
	

	//Save all accumulated star data
	//Save results to file
	ContructFilePath("Stars.txt", "Output", filename); //file where data about stars should go
	if(NStarsCurrent>0){
		SaveStarMatrix(filename, IsEnd);
		FreeStarMatrix(); //free memory
	}
	//Save accumulated core data
	ContructFilePath("Clouds.txt", "Output", filename); //file where data about clouds should go
	for ( l = 0; l < CMFTIMENUM; l++ ){ //check all CMF times 
		CMFFileName(Corefilename, filename, CMFTimes[l]);
		if(NCoresCurrent[l]>0){
			SaveCoreMatrix(l,Corefilename,IsEnd);
		}
	}
	FreeCoreMatrix(); //free memory
	
}


//Check if memory limit for star and core matrices is exceeded
bool TooLargeStarsAndCoresMatrices(){
	long TotalAlloc=NStarsAllocated*8; //stars, 8 data items per star
	long l;
	
	for ( l = 0; l < CMFTIMENUM; l++ ){ //check all CMF times
		TotalAlloc+=NCoresAllocated[l]*12; //12 numbers per core
	}	
	if(TotalAlloc<MaxAllocated){
		return false;
	}
	else
	{
		return true;
	}	
}

/*************************************************************/
		//Saving to disk
		
		
//Save StarMatrix
void SaveStarMatrix(char* filename, bool IsEnd){
	long i,j;
	ofstream myfile; //Output file handle
    //Open file
    myfile.open(filename,ios::app);
    //Save parameters of stars
    for(i=0;i<NStarsCurrent;i++){ 
		for(j=0;j<8;j++){
			if(j<7){
				myfile<<StarMatrix[j][i]<<"\t";
			}
			else{
				myfile<<StarMatrix[j][i]<<"\n";
			}
     	}
     }
     if(IsEnd){
     	myfile <<"END"<<"\n"; //mark end for this parent cloud in file
	 }
     myfile.flush(); //finish writing
     myfile.close(); //close file	
}

//Save CoreMatrix
void SaveCoreMatrix(long timeindex, char* filename, bool IsEnd){
	long i,j;
	ofstream myfile; //Output file handle
    //Open file
    myfile.open(filename,ios::app);
    //Save parameters of cores for this obs. time
    for(i=0;i<NCoresCurrent[timeindex];i++){ 
		for(j=0;j<12;j++){
			if(j<11){
				myfile<<CoreMatrix[timeindex][j][i]<<"\t";
			}
			else{
				myfile<<CoreMatrix[timeindex][j][i]<<"\n";
			}
     	}
     }
     if(IsEnd){
     	myfile <<"END"<<"\n"; //mark end for this parent cloud in file
	 }
     myfile.flush(); //finish writing
     myfile.close(); //close file	
}


//Save field data
void SaveFields(struct densitystruct_s* fieldstruct, long ObjectTimeindex){
	long i;
	char filename[200],protofilename[200];
	ofstream myfile; //Output file handle
	
	sprintf(protofilename, "fields_t%g.txt", CMFTimes[ObjectTimeindex]);
	ContructFilePath(protofilename, "Output", filename); //file where data about fields should go
	
    //Open file
    myfile.open(filename,ios::app);
	#if (ExtraVerboseOutput==true)
	cout<<"Fields output file "<< filename<<"\n";
	#endif
    //Save field data
    for ( i = 0; i < NgridTotal; i++ ){
		myfile<<fieldstruct[0].xcoord[i]<<"\t"<<fieldstruct[0].ycoord[i]<<"\t"<<fieldstruct[0].zcoord[i]<<"\t"<<fieldstruct[0].rho_field[i]<<"\t"<<fieldstruct[0].T_field[i]<<"\t"<<fieldstruct[0].volume<<"\n";
     }
     myfile.flush(); //finish writing
     myfile.close(); //close file	
}

/*************************************************************/
		//Storing data
		
//Store parameters of a protostar		
void AssignStarMatrixValues(double Mass, double X0, double Y0, double Z0, double RelativeVelocityX, double RelativeVelocityY, double RelativeVelocityZ, double tcurrent){
	StarMatrix[0][NStarsCurrent]=Mass;
	StarMatrix[1][NStarsCurrent]=X0;
	StarMatrix[2][NStarsCurrent]=Y0;
	StarMatrix[3][NStarsCurrent]=Z0;
	StarMatrix[4][NStarsCurrent]=RelativeVelocityX;
	StarMatrix[5][NStarsCurrent]=RelativeVelocityY;
	StarMatrix[6][NStarsCurrent]=RelativeVelocityZ;
	StarMatrix[7][NStarsCurrent]=tcurrent;
	NStarsCurrent++; //one more star has been added
}

//Store parameters of a core		
void AssignCoreMatrixValues(long timeindex, double MassInit, double Mass, double X0, double Y0, double Z0, double RelativeVelocityX, double RelativeVelocityY, double RelativeVelocityZ, double Machedge, double Temperature, double Radius, double RhoAverage){
	CoreMatrix[timeindex][0][NCoresCurrent[timeindex]]=MassInit;
	CoreMatrix[timeindex][1][NCoresCurrent[timeindex]]=Mass;
	CoreMatrix[timeindex][2][NCoresCurrent[timeindex]]=X0;
	CoreMatrix[timeindex][3][NCoresCurrent[timeindex]]=Y0;
	CoreMatrix[timeindex][4][NCoresCurrent[timeindex]]=Z0;
	CoreMatrix[timeindex][5][NCoresCurrent[timeindex]]=RelativeVelocityX;
	CoreMatrix[timeindex][6][NCoresCurrent[timeindex]]=RelativeVelocityY;
	CoreMatrix[timeindex][7][NCoresCurrent[timeindex]]=RelativeVelocityZ;
	CoreMatrix[timeindex][8][NCoresCurrent[timeindex]]=Machedge;
	CoreMatrix[timeindex][9][NCoresCurrent[timeindex]]=Temperature;
	CoreMatrix[timeindex][10][NCoresCurrent[timeindex]]=Radius;
	CoreMatrix[timeindex][11][NCoresCurrent[timeindex]]=RhoAverage;
	NCoresCurrent[timeindex]++; //one more star has been added
}


/**************************************************************/
		//Init
		
void InitStarsAndCores()
{
	
	long i;
	
	//Init global counting variables
	NStarsAllocated=0;
	NStarsCurrent=0;
	NCoresAllocated=(long*)calloc(CMFTIMENUM, sizeof (long));
	NCoresCurrent=(long*)calloc(CMFTIMENUM, sizeof (long));
	StarMatrix=(double**)calloc(8,sizeof(double*));
	CoreMatrix=(double***)calloc(CMFTIMENUM,sizeof(double**));
	for(i=0;i<CMFTIMENUM;i++){
		CoreMatrix[i]=(double**)calloc(12,sizeof(double*));
	}

}

// Constructs file name for CMF evolution file
void CMFFileName(char* fullname, char* filename, double CMFTime)
{
	char ending[100];
	char beginning[100];
	char* chpointer;

	//find file extension at the end of the string
	chpointer=strchr(filename,'x');
	//make a copy of the filename without extension
	strncpy(beginning, filename, chpointer-filename-2);
	beginning[chpointer-filename-1]=0;
	sprintf( ending, "_cores_t%g.txt", CMFTime);
	fullname[0]=0; //set it as empty string
	strcat (fullname,beginning);
	strcat (fullname,ending);	
}

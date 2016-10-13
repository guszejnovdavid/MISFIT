/**************************************************************************************/
			/*  Handling Input data structures for Evolvefield */
			
// Collection of routines used by Master for storing input data

#include <stdlib.h>
#include <stddef.h>
#include "vectorfunctions.h"
#define InputDataManage_DEFINED // tells StarsAndCores.h not to define the external variables as we will do it here
#include "GlobalVars.h"
#include "InputDataManage.h"
#include <fstream>
#include <iostream>
#include "EvolveField.h"
# include "mpi.h"
#include "FileManage.h"

/*****************************************/
		//Global variables
		
MPI_Datatype mpi_inputstruct_type; //declare the new mpi type of the structure
MPI_Datatype mpi_densitystruct_type; //declare the new mpi type of the density structure  		
long NDataPerStruct=2*NgridTotal+10;
struct  inputstruct_s  *InputStorage; //vector, containing the input parameter stuctures used for calling new runs.
double* InputPriority; //vector containing the priority of each run
bool* IsUsed; // vector showing which allocated spaces are used
long NInputAllocated=0; //Number of structures allocated in InputStorage
long NInputCurrent=0; //Number of structures stored in InputStorage


/*************************************************************/
		//Changing size
		
		
//Increases the number of structs by StructIncreaseUnit
void IncreaseInputStorage(){
	long i;
	InputStorage=(struct inputstruct_s*) realloc(InputStorage, (NInputAllocated+StructIncreaseUnit)* sizeof(struct inputstruct_s));
	if(InputStorage==NULL){ //unsuccessful allocation
		cout<<"Memory allocation failure! Exiting!";
		exit(EXIT_FAILURE);
	}
	IsUsed=(bool*) realloc(IsUsed, (NInputAllocated+StructIncreaseUnit)*sizeof(bool));
	for(i=0;i<StructIncreaseUnit;i++){
		IsUsed[NInputAllocated+i]=false;
	}
    InputPriority=(double*) realloc(InputPriority, (NInputAllocated+StructIncreaseUnit)*sizeof(double));
	NInputAllocated+=StructIncreaseUnit; //note the number of elements
	cout<<"Input storage expanded to "<<NInputAllocated<<"\n";
}

//Free StarMatrix
void FreeInputStorage(){
	NInputAllocated=0;
	NInputCurrent=0;
	free(InputStorage);
	free(IsUsed);
	free(InputPriority);
}


//Save input data and free memory
void SaveAndFreeInputStorage(){
	char filename[200];
	
	//Get filename
	ContructFilePath("Restart_Inputs.txt", "Output", filename);
	//Save all accumulated input data
	SaveInputStorage(filename);
	FreeInputStorage(); //free memory
}


/*************************************************************/
		//Saving to disk
		
		
//Save StarMatrix
void SaveInputStorage(char* filename){
	long i,j;
	ofstream myfile; //Output file handle
	
	if(file_exists (filename)){
		remove(filename); //removes previous restart data
	}
    //Open file
    myfile.open(filename,ios::out);
    //Save parameters of stars
    for(i=0;i<NInputAllocated;i++){
    	if(IsUsed[i]){ //only save usable data
	    	//Save rho
			for(j=0;j<NgridTotal;j++){
				myfile<<InputStorage[i].rho_field[j]<<"\n";
	     	}
	    	//Save T
			for(j=0;j<NgridTotal;j++){
				myfile<<InputStorage[i].T_field[j]<<"\n";
	     	}
	     	//Save R0
	     	myfile<<InputStorage[i].R0<<"\n";
	     	//Save X0,Y0,Z0
	     	myfile<<InputStorage[i].X0<<"\n"; myfile<<InputStorage[i].Y0<<"\n";	myfile<<InputStorage[i].Z0<<"\n";
	     	//Save Machedge
	     	myfile<<InputStorage[i].Machedge<<"\n";
	     	//Save VX0,VY0,VZ0
	     	myfile<<InputStorage[i].RelativeVelocityX<<"\n"; myfile<<InputStorage[i].RelativeVelocityY<<"\n";	myfile<<InputStorage[i].RelativeVelocityZ<<"\n";
	     	//Save InitScale
	     	myfile<<InputStorage[i].InitScale<<"\n";
	     	//Save tparent
	     	myfile<<InputStorage[i].tparent<<"\n";  // new line after last data item
		}
     }
     myfile.flush(); //finish writing
     myfile.close(); //close file	
}


/*************************************************************/
		//Storing data
		
//Store parameters of a protostar		
void PutIntoInputStorage(struct inputstruct_s* inputstruct, long pos_index){
	long i;
	
	//Rho
	for(i=0;i<NgridTotal;i++){
		InputStorage[pos_index].rho_field[i]=inputstruct[0].rho_field[i];
    }
	//T
	for(i=0;i<NgridTotal;i++){
		InputStorage[pos_index].T_field[i]=inputstruct[0].T_field[i];
    }
    InputStorage[pos_index].R0=inputstruct[0].R0;
    InputStorage[pos_index].X0=inputstruct[0].X0;
    InputStorage[pos_index].Y0=inputstruct[0].Y0;
    InputStorage[pos_index].Z0=inputstruct[0].Z0;
    InputStorage[pos_index].Machedge=inputstruct[0].Machedge;
    InputStorage[pos_index].RelativeVelocityX=inputstruct[0].RelativeVelocityX;
    InputStorage[pos_index].RelativeVelocityY=inputstruct[0].RelativeVelocityY;
    InputStorage[pos_index].RelativeVelocityZ=inputstruct[0].RelativeVelocityZ;
    InputStorage[pos_index].InitScale=inputstruct[0].InitScale;
    InputStorage[pos_index].tparent=inputstruct[0].tparent;
    
    //Assigned values
    IsUsed[pos_index]=true; //space occupied
    InputPriority[pos_index]=PriorityFunction(inputstruct);
    //InputPriority[pos_index]=1.0/inputstruct.R0;
    
	NInputCurrent++; //one more input struct has been added
}


/**************************************************************/
		//Initialize input data type and make it MPI passable
		


void InitInputDataType(){

//MPI Variables
const int nitems=12; //number of items in the structure	
int blocklengths[12]={NgridTotal,NgridTotal,1,1,1,1,1,1,1,1,1,1}; //length of struct blocks
MPI_Datatype types[12] = {MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE}; //type of variables in blocks
MPI_Aint offsets[12]; //memory offsets for each block	

//Find offsets of input datatype
offsets[0] = offsetof(inputstruct_s,rho_field);
offsets[1] = offsetof(inputstruct_s, T_field);
offsets[2] = offsetof(inputstruct_s, R0);
offsets[3] = offsetof(inputstruct_s, X0);
offsets[4] = offsetof(inputstruct_s, Y0);
offsets[5] = offsetof(inputstruct_s, Z0);
offsets[6] = offsetof(inputstruct_s, Machedge);
offsets[7] = offsetof(inputstruct_s, RelativeVelocityX);
offsets[8] = offsetof(inputstruct_s, RelativeVelocityY);
offsets[9] = offsetof(inputstruct_s, RelativeVelocityZ);
offsets[10] = offsetof(inputstruct_s, InitScale);
offsets[11] = offsetof(inputstruct_s, tparent);


//create MPI passable struct type
MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_inputstruct_type); 
//allows the type to be used for MPI communication
MPI_Type_commit(&mpi_inputstruct_type); 

	
}

/**************************************************************/
		//Initialize density structure data type and make it MPI passable
		


void InitDensityStructDataType(){
//MPI Variables
const int nitems=6; //number of items in the structure	
int blocklengths[6]={NgridTotal,NgridTotal,NgridTotal,NgridTotal,NgridTotal,1}; //length of struct blocks
MPI_Datatype types[6] = {MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE}; //type of variables in blocks
MPI_Aint offsets[6]; //memory offsets for each block	

//Find offsets of input datatype
offsets[0] = offsetof(densitystruct_s,rho_field);
offsets[1] = offsetof(densitystruct_s, T_field);
offsets[2] = offsetof(densitystruct_s, xcoord);
offsets[3] = offsetof(densitystruct_s, ycoord);
offsets[4] = offsetof(densitystruct_s, zcoord);
offsets[5] = offsetof(densitystruct_s, volume);



//create MPI passable struct type
MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_densitystruct_type); 
//allows the type to be used for MPI communication
MPI_Type_commit(&mpi_densitystruct_type); 

	
}

/**************************************************************/
		//Initialize storage of input srtuctures

void InitInputDataStorage(long Nitem)
{
	//Init global counting variables
	NInputAllocated=Nitem;
	NInputCurrent=0;
	InputStorage=(struct inputstruct_s*) calloc(Nitem,sizeof(struct inputstruct_s));
	if(InputStorage==NULL){ //unsuccessful allocation
		cout<<"Memory allocation failure! Exiting!";
		exit(EXIT_FAILURE);
	}
    IsUsed=vector_bool(Nitem);
    InputPriority=vector_double(Nitem);
}

/***************************************************/
		//Calculate priority for input dataset
double PriorityFunction(struct inputstruct_s* inputstruct){
	//This function sets the order in which the input datasets are being handed out to the slaves. To conserve memory, it is better to call runs that will
	//not create new data structs (e.g. terminate without fragmenting)
	double res;
	
	res=1.0/inputstruct[0].R0;
	
	return res;
}


/***************************************************/
//Finds the index of the dataset with the highest priority
long NextToCall(){
	long i,ind=-1;
	double temp=0;
	
	for(i=0;i<NInputAllocated;i++){
		if(IsUsed[i]){
			if(ind==-1){
				ind=i;
				temp=InputPriority[i];
			}
			else{
				if(InputPriority[i]>temp){
					ind=i;
					temp=InputPriority[i];
				}
			}
			
		}
    }
    
    return ind;	
}

/**************************************************/
//Finds first unused slot in the dataset
long FirstUnused(){
	long i;
	for(i=0;i<NInputAllocated;i++){
		if(IsUsed[i]==false){
			return i;
		}
	}
	cout<<"All of InputStorage used.";
	return -1; //error
}

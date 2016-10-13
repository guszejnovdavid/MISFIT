/**************************************************************************************/
								/*  MISFIT Restart Routines */
								
// Routines needed to save and restart calculation

#include "Restart.h"
#include "vectorfunctions.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "FileManage.h"
#include "GlobalVars.h"
#include "InputDataManage.h"

using namespace std;


//Read restart file data
void ReadRestartFile(long NData)
{

	long i,j;
	char line[100];
	FILE* sourcefile;
    char filename[200];

	//Get filename
	ContructFilePath("Restart_Inputs.txt", "Output", filename);
	//open file
    sourcefile=fopen(filename,"r");
    for(i=0;i<NData;i++){
    	// Rho
    	for(j=0;j<NgridTotal;j++){
     		fgets (line, 100, sourcefile);
			InputStorage[i].rho_field[j]=strtod(line, NULL);
		}
    	// T
    	for(j=0;j<NgridTotal;j++){
     		fgets (line, 100, sourcefile);
			InputStorage[i].T_field[j]=strtod(line, NULL);
		}
		// R0
     	fgets (line, 100, sourcefile);
		InputStorage[i].R0=strtod(line, NULL);
		// X0,Y0,Z0
     	fgets (line, 100, sourcefile);
		InputStorage[i].X0=strtod(line, NULL);
     	fgets (line, 100, sourcefile);
		InputStorage[i].Y0=strtod(line, NULL);
     	fgets (line, 100, sourcefile);
		InputStorage[i].Z0=strtod(line, NULL);
		// MachEdge
     	fgets (line, 100, sourcefile);
		InputStorage[i].Machedge=strtod(line, NULL);
		// VX0,VY0,VZ0
     	fgets (line, 100, sourcefile);
		InputStorage[i].RelativeVelocityX=strtod(line, NULL);
     	fgets (line, 100, sourcefile);
		InputStorage[i].RelativeVelocityY=strtod(line, NULL);
     	fgets (line, 100, sourcefile);
		InputStorage[i].RelativeVelocityZ=strtod(line, NULL);
		// Initscale
     	fgets (line, 100, sourcefile);
		InputStorage[i].InitScale=strtod(line, NULL);
		// tparent
     	fgets (line, 100, sourcefile);
		InputStorage[i].tparent=strtod(line, NULL);	
		//Vector elements
		IsUsed[i]=true;
		InputPriority[i]=PriorityFunction(&InputStorage[i]);
    }
	fclose(sourcefile);
	NInputCurrent=NData;
	remove(filename);
		
}


//Saves loop control data
void SaveRestartControl(long Niter){
	ofstream myfile; //Output file handle
	char filename[200];
	
	//Get filename
	ContructFilePath("RestartControl.txt", "Output", filename);
	//open output file
	myfile.open(filename,ios::out);
    //Save parameters
    myfile<<Niter<<"\n";
    myfile.flush(); //finish writing
    myfile.close(); //close file
}

//Saves loop control data
void ReadRestartControl(long* Niter){
	char filename[200];
	FILE *sourcefile; //Input file handle
	char line[100];
	
	//Get filename
	ContructFilePath("RestartControl.txt", "Output", filename);
	//open input file
	sourcefile=fopen (filename,"r");
    //Read parameters
    fgets (line, 100, sourcefile);
    Niter[0]=strtol(line,NULL, 0);
    fclose(sourcefile);
    remove(filename); //remove file
}


// Collection of routines to find files, create directories (cross platform compatible)

#include <stdio.h>
#include "FileManage.h"

#ifdef __linux__
  // Linux Includes Here
  #include <sys/types.h>
  #include <sys/stat.h>
  const char SlashOS='/';
#endif

#if defined(_WIN32) || defined(_WIN64)
  // Windows Includes Here
  #include <windows.h>
  const char SlashOS='\\';
#endif

//Checks if a file exists
bool file_exists (char* filename) {
    if (FILE *file = fopen(filename, "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}

//Construct filepath for file in a directory
void ContructFilePath(char* filename, char* directoryname, char* fullpath){	
	sprintf( fullpath, "%s%c%s", directoryname,SlashOS,filename );
}


/**********************************************************/
		/* Create directory if it doesn't exist */
		
//WINDOWS
#if defined(_WIN32) || defined(_WIN64)		
	int CreateDir(char* directoryname){
	   return CreateDirectory (directoryname, NULL); //returns 0 if successful
	}
#endif

//LINUX
#ifdef __linux__
	int CreateDir(char* directoryname){
	   return mkdir(directoryname, 0770);; //returns 0 if successful
	}
#endif

/***************************************************************/
	/*  Writing a line to a text to file */
	
	
void WriteLineToTextFile(char* line, char* filename)
{
	FILE *targetfile; //input file handle
    //Open file
    targetfile = fopen (filename,"a");
    //Save line to file
	fprintf(targetfile,"%s\n",line);
	fflush(targetfile);
	fclose(targetfile);		
}


//Count the number of lines in a file
long CountNumberOfLines(char* filename){
	
FILE *fp = fopen(filename,"r");
int ch=0;
long lines=0;

while(!feof(fp))
{
  ch = fgetc(fp);
  if(ch == '\n')
  {
    lines++;
  }
}

return lines;
}


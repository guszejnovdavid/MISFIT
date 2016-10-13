//Checks if a file exists
bool file_exists (char* filename);

//Construct filepath for file in a directory
void ContructFilePath(char* filename, char* directoryname, char* fullpath);

/**********************************************************/
		/* Create directory if it doesn't exist */	
int CreateDir(char* directoryname);

/***************************************************************/
	/*  Writing a line to a text to file */
void WriteLineToTextFile(char* line, char* filename);

//Count the number of lines in a file
long CountNumberOfLines(char* filename);



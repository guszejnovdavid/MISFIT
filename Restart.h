#ifndef Restart_INCLUDED
#define Restart_INCLUDED

void ReadRestartFile(long NData);
void SaveRestartControl(long Niter);
void ReadRestartControl(long* Niter);

#endif

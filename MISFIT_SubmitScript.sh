### PBS submit script for MISFIT
#!/bin/sh
### Job name job (default is name of pbs script file)
#PBS -N M1E3_dens_isoT
### Number and type of nodes:
#PBS -l nodes=10:core12
### Time limit([[h:]m:]s).
#PBS -l walltime=48:00:00
#PBS -mea
### Join output and error files 
#PBS -j oe
### Output file
#PBS -o MISFIT.out
#PBS -V
### Email me at this address.
#PBS -M guszejnov@caltech.edu
### Send me email when job begins (b), aborts (a) or ends (e)
#PBS -m abe

Ntries=1

cd /panfs/ds08/hopkins/guszejnov/MISFIT/density_test/M1E3_dens_isoT
mpirun -np 120 MISFIT $Ntries
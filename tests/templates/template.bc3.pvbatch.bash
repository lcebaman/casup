#!/bin/bash
#$Id: template.bc3.pvbatch.bash 175 2015-12-15 12:31:30Z mexas $
#PBS -l walltime=01:00:00,nodes=1:ppn=16
#PBS -j oe
#PBS -q testq

cd $HOME/nobkp/cgpack/head/tests
module load apps/paraview-4.3.1
mpirun -n 16 xvfb-run pvbatch --use-offscreen-rendering zcracks.py 
#xvfb-run pvbatch --use-offscreen-rendering zgb.py 
#mpirun -n 16 pvbatch --use-offscreen-rendering z.xdmf.py
#mpirun -n 16 pvbatch zgb.py
#mpirun -n 16 pvbatch z.xdmf.py

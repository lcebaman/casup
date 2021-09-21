#!/bin/bash --login
#PBS -l mppwidth=64
#PBS -l mppnppn=32
#PBS -l walltime=00:05:00
#PBS -A e277
#PBS -j oe

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
export WORK=`echo $HOME | sed s/home/work/g`
export XDG_CONFIG_HOME=$WORK/ParaViewIniConfig

module load paraview-servers/3.12.0

if [ -z ${PARAVIEW_SERVER_DIR} ] ; then
  echo "Error: PARAVIEW_SERVER_DIR not set. Exiting"
  exit 4
fi

echo "DEBUG: PBS_O_WORKDIR is" $PBS_O_WORKDIR
echo "DEBUG: PARAVIEW_SERVER_DIR is" ${PARAVIEW_SERVER_DIR}
echo "DEBUG: XDG_CONFIG_HOME is" ${XDG_CONFIG_HOME}

export MPPWIDTH=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
export MPPNPPN=`qstat -f $PBS_JOBID | awk '/mppnppn/ {print $3}'`

aprun -n ${MPPWIDTH} -N ${MPPNPPN} ${PARAVIEW_SERVER_DIR}/bin/pvbatch \
       --use-offscreen-rendering \
       pv.py

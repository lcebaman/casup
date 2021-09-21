#!/bin/bash --login
#
#$Id: hec.pvserver.sample 191 2015-12-15 21:46:16Z mexas $
#
#PBS -N pvserver
#PBS -l mppwidth=32
#PBS -l mppnppn=32
#PBS -l walltime=01:00:00
#PBS -A e277
#PBS -v HOST_IP
#PBS -m a

export XDG_CONFIG_HOME=${WORK}/ParaViewIniConfig 
cd $PBS_O_WORKDIR
# this job script is to launch the paraview server as a parallel job
# on one node with 32 cores. It will connect with client running on HOST_IP

module load paraview-servers/3.12.0

if [ -z ${PARAVIEW_SERVER_DIR} ] ; then
  echo "Error: PARAVIEW_SERVER_DIR not set. Exiting"
  exit 4
fi

echo "DEBUG: PARAVIEW_SERVER_DIR is" ${PARAVIEW_SERVER_DIR}
echo "DEBUG: expect to connect to client on host ip address= ${HOST_IP}"

export MPPWIDTH=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
export MPPNPPN=`qstat -f $PBS_JOBID | awk '/mppnppn/ {print $3}'`

aprun -n ${MPPWIDTH} -N ${MPPNPPN} ${PARAVIEW_SERVER_DIR}/bin/pvserver \
       --use-offscreen-rendering \
       --reverse-connection \
       --server-port=75000  \
       --client-host=${HOST_IP}

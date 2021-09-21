#!/bin/bash
# submit from ./bin.scorep directory with "sbatch scorep.sbatch"
#
#SBATCH --job-name="CGPACK"	# Job name
#SBATCH --output=CGPACK-%A.out	# Job output file
#SBATCH --ntasks=32		# Total number of MPI processes
#SBATCH --cpus-per-task=1	# Number of cores for each MPI process (OpenMP threads)
#SBATCH --time=00:05:00		# Maximum run time (hh:mm:ss)
#SBATCH --exclusive		# No sharing of compute nodes
#SBATCH --partition=sandybridge # Darwin compute nodes
#SBATCH --account=TRAINING
#SBATCH --reservation=PAW	# Dedicated workshop reservation

# set up environment
. /etc/profile.d/modules.sh    
module purge
module load default-impi-LATEST

# benchmark configuration
export NPB_MZ_BLOAD=0
export OMP_NUM_THREADS=4
PROCS=32
CLASS=B
EXE=./bt-mz_$CLASS.$PROCS

# Score-P measurement configuration
#export SCOREP_EXPERIMENT_DIRECTORY=scorep_bt-mz_sum
#export SCOREP_FILTERING_FILE=../config/scorep.filt
#export SCOREP_METRIC_PAPI=PAPI_TOT_INS,PAPI_FP_INS
#export SCOREP_METRIC_PAPI_PER_PROCESS=PAPI_L2_TCM
#export SCOREP_METRIC_RUSAGE=ru_stime
#export SCOREP_METRIC_RUSAGE_PER_PROCESS=ru_maxrss
#export SCOREP_TIMER=gettimeofday
#export SCOREP_TOTAL_MEMORY=21MB

set -x
EXECUTABLE=testABW.x
CACONFILE=xx14.conf

# Prepare conf file
echo -genvall -genv I_MPI_FABRICS=shm:dapl \
        -n $PROCS ./$EXECUTABLE > $CACONFILE

$EXECUTABLE

# -machinefile ./$confile \

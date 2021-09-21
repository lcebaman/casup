#!/bin/sh

#$Id: instrument.sh 101 2016-05-17 15:30:37Z mexas $
#
# Copyright (c) 2016 Anton Shterenlikht, The University of Bristol
#
# See LICENSE for licensing conditions

export TAU_OPTIONS="-optShared -optVerbose -optCompInst"
#export TAU_MAKEFILE=$HOME/tau-2.25.1.1/x86_64/lib/Makefile.tau-icpc-papi-mpi-pdt
export TAU_MAKEFILE=$HOME/tau-2.25.2-intel/x86_64/lib/Makefile.tau-icpc-papi-mpi-pdt

make clean -i -f Makefile-bc3-mpiifort-tau
make all deinstall install -i -f Makefile-bc3-mpiifort-tau

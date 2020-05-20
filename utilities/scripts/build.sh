#!/bin/bash

set -v

#uname -s
#uname -n
#
which gfortran
which gcc
which g++

#make RFLU=1

if [ -f "openmpi-2.0.1/config.log" ]; then
	echo "Using cached OpenMPI"
	echo "Configuring OpenMPI"
	cd openmpi-2.0.1
	./configure --prefix=$TRAVIS_BUILD_DIR/openmpi CC=gcc CXX=g++ FC=gfortran &> openmpi.configure
else
# install OpenMPI from source
	echo "Downloading OpenMPI Source"
	wget https://www.open-mpi.org/software/ompi/v2.0/downloads/openmpi-2.0.1.tar.gz
	tar zxf openmpi-2.0.1.tar.gz
	echo "Configuring and building OpenMPI"
	cd openmpi-2.0.1
	echo "In dir"
	./configure --prefix=$TRAVIS_BUILD_DIR/openmpi CC=gcc CXX=g++ FC=gfortran &> openmpi.configure
	echo "After configure"
	make -j4 &> openmpi.make
	make install &> openmpi.install
	cd ..
fi

# recommended by Travis CI documentation to unset these for MPI builds
test -n $CC && unset CC
test -n $CXX && unset CXX
test -n $FC && unset FC

which mpif90
which mpic++

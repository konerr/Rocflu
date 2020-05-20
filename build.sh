#!/bin/bash

set -ev

uname -s
uname -n

#ls $TRAVIS_BUILD_DIR/openmpi/bin/mpi*
#echo "After ls"
#sudo cp $TRAVIS_BUILD_DIR/openmpi/bin/mpi* /usr/bin/
#echo "After cp"

ls /usr/
#ls /usr/lib/*metis*
find /usr/ -name "*metis*" 

which gfortran
which cpp
#which mpif90
#which mpic++

#make clean
#echo 'Cleaning directory'
#echo '--------------------------'
#echo 'Runnnig RFLU=1'
#make RFLU=1


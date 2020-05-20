#!/bin/bash

set -ev

#uname -s
#uname -n
#
which mpif90
which mpic++
which gfortran
which cpp

make clean
echo 'Cleaning directory'
echo '--------------------------'
echo 'Runnnig RFLU=1'
make RFLU=1


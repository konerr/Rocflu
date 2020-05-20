#!/bin/bash

set -v

uname -s
uname -n
#
which gfortran
which cpp
which mpif90
which mpic++

make clean
echo 'Cleaning directory'
echo '--------------------------'
echo 'Runnnig RFLU=1'
make RFLU=1


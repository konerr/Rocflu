!#/bin/bash

set -ev

uname -s
uname -n

which mpif90
which mpic++
which gfortran
which cpp

make RFLU=1

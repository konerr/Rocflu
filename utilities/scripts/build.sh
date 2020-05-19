!#/bin/bash

set -ev

uname -s
uname -n

which mpif90
which mpic++

make RFLU=1

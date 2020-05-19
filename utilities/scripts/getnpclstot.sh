#!/bin/bash

################################################################################
#
# Purpose: Script for extracting the total number of particles in a given 
#   RocfluMP computation from the .pdim files of each region. 
#
# Description: The script determines the number of regions (and hence the 
#   number of pdim files) from the .map file and then proceeds to read the 
#   .pdim files of each region for the given time stamp to compute the total
#   number of particles across all regions for that time stamp. 
#
# Notes: The time stamp provided to this script is assumed to be in the right
#   format so it can be used immediately - i.e., it is identical to the time
#   stamp in the .pdim file name.
# 
# Author: Andreas Haselbacher
#
################################################################################

if [ $# -lt 2 ]; then
  echo "Usage: $0 <casename> <time stamp>"
  exit 1
fi

# Set mapping file name and check whether it exists
fileNameMap=$1.map

if [ -f $fileNameMap ]; then
  exec 3< $fileNameMap
else 
  echo $fileNameMap "file does not exist."
  exit 1
fi

# Read header lines in map file and find nRegions
read <&3 line
read <&3 line
read <&3 nRegions

# Initialize loop variables
iReg=0
nPclsTot=0

# Loop over regions, build pdim file name, read nPcls, and accumulate
while [ $iReg -lt $nRegions ]; do 
  let iReg=iReg+1   
  
  if [ $iReg -lt 10 ]; then
    fileNameDim=$1.pdim_0000$iReg\_$2
  elif [ $iReg -lt 100 ]; then
    fileNameDim=$1.pdim_000$iReg\_$2
  elif [ $iReg -lt 1000 ]; then
    fileNameDim=$1.pdim_00$iReg\_$2
  elif [ $iReg -lt 10000 ]; then
    fileNameDim=$1.pdim_0$iReg\_$2
  else
    echo "ERROR - Too many files for this script, it needs to be modified."
    exit 1
  fi

  if [ ! -f $fileNameDim ]; then
    echo $fileNameDim "file does not exist."
    exit 1
  fi 

  exec 3< $fileNameDim
 
  read <&3 line
  read <&3 line
  read <&3 nPcls

  let nPclsTot=nPclsTot+$nPcls
done

echo $nPclsTot

exit 0

################################################################################
#
# RCS Revision history:
#
# $Log: getnpclstot.sh,v $
# Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
# merged rocflu micro and macro
#
# Revision 1.1.1.1  2014/07/15 14:31:37  brollin
# New Stable version
#
# Revision 1.1  2007/04/12 19:02:51  haselbac
# Initial revision
#
################################################################################

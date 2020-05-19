#!/bin/bash

################################################################################
#
# Purpose: Wrapper script for getnpclstot.sh script. 
#
# Description: The script processes the file <filename> line by line and calls
#   getnpclstot.sh once for each line. The lines in file <filename> are assumed 
#   to consist only of numeric strings representing time stamps at which the 
#   the solution is stored. The command cut is used to strip out the part of 
#   the time stamps which is actually encoded into file names.
#
# Notes: Based on script linebyline.sh found on web.
# 
# Author: Andreas Haselbacher
#
################################################################################

if [ $# -lt 1 ]; then
  echo "Usage: $0 <casename> < <filename>"
  exit 1
fi

while true
do
  read line

  if [ $? -eq 0 ]; then

# Get position of "E" in and length of time string
    ep=$(echo $line | sed -n "s/[E].*//p" | wc -c)  
    ll=${#line}
 
    epp1=$[$ep + 1]
    epm1=$[$ep - 1]
    llm1=$[$ll - 1]

# Strip out various parts of time string
    timeStamp1=$(echo $line | cut -c1-${epm1})
    timeStamp2=$(echo $line | cut -c${ep}-${epp1})
    timeStamp3=$(echo $line | cut -c${llm1}-${ll})

# Round time stamp (without exponent) to 5 digits
    timeStamp1r=$(printf "%1.5f\n" $timeStamp1)

# Put string back together
    timeStamp=${timeStamp1r}${timeStamp2}${timeStamp3}

    echo $line $(getnpclstot.sh $1 $timeStamp)
  else
    exit 1
  fi
done

################################################################################
#
# RCS Revision history:
#
# $Log: getnpclstotwrap.sh,v $
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

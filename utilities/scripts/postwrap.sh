#!/bin/bash

################################################################################
#
# Purpose: Wrapper script for rflupost code. 
#
# Description: The script processes the file <filename> line by line and calls
#   rflupost once for each line. The lines in file <filename> are assumed to 
#   to consist only of numeric strings. 
#
# Notes: Based on script linebyline.sh found on web.
# 
# Author: Andreas Haselbacher
#
################################################################################

if [ $# -lt 2 ]; then 
  echo "Usage: postwrap <casename> <verbosity> < <filename>"
  exit 1
fi

while true
do
  read line
  
  if [ $? -eq 0 ]; then
    /home/brollin/rocfluMacro_11Aug2014_postDebug/rflupost -c $1 -v $2 -s $line
  else
    exit 0
  fi
done

################################################################################
#
# RCS Revision history:
#
# $Log: postwrap.sh,v $
# Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
# merged rocflu micro and macro
#
# Revision 1.1.1.1  2014/07/15 14:31:37  brollin
# New Stable version
#
# Revision 1.1  2007/04/12 19:02:42  haselbac
# Initial revision after split from RocfloMP
#
# Revision 1.2  2005/08/23 16:36:10  haselbac
# Adapted to new rflupost
#
# Revision 1.1  2004/10/20 12:58:24  haselbac
# Initial revision
#
################################################################################

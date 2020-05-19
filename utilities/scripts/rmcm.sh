#!/bin/sh

################################################################################
#
# Purpose: Remove carriage return (ctrl M) from files. 
#
# Description: None.
#
# Notes: 
#   1. IMPORTANT: This file may not be displayed properly in CVSWeb because the
#      ctrl-M character is interpreted. Nevertheless, the file is stored 
#      correctly in CVS, so DO NOT CHANGE IT.
# 
# Author: Andreas Haselbacher
#
################################################################################

TEMPFILE=rmcm.tmp

sed 's///g' $1 > $TEMPFILE
mv $TEMPFILE $1

################################################################################
#
# RCS Revision history:
#
# $Log: rmcm.sh,v $
# Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
# merged rocflu micro and macro
#
# Revision 1.1.1.1  2014/07/15 14:31:37  brollin
# New Stable version
#
# Revision 1.1  2007/04/12 19:02:42  haselbac
# Initial revision after split from RocfloMP
#
# Revision 1.2  2005/09/19 18:27:38  haselbac
# Added note
#
# Revision 1.1  2005/09/19 16:30:48  haselbac
# Initial revision
#
################################################################################

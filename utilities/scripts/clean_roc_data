#!/bin/bash

################################################################################
#
# Purpose: To clean up data from a Rocflu run to start from scratch.
#
#
# Description: This code determines how many solution files are present in the
#              directory and runs the Rocflu post processing routine,rflupost,
#              on the files in the directory. It also creates a file called
#              "Extracted_Data_Labels.txt" that holds all of the timestamps
#              of the solution files that were processed.
#
# Example: clean_rocflu_data <casename> 
#	   <casename>= Identifying casename that's on Rocflu data i.e. <casename>.map
#
#
#
#
# Author: Christopher Neal
# Date:   05/22/2014
# Updated: 03/27/2015
################################################################################


rm -rf $1.com_* $1.fom_* $1.cmp* $1.mixt* $1.pcoa* $1.dim* $1.map $1.rnm* $1.rin $1.grda* rflu* $1.con *.log $1.ptm $1.pm

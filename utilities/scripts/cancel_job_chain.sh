#!/bin/bash
#######################################################################################
# The purpose of this script is to cancel several depdendent jobs. This file looks at
# the file which contains all of the jobIDs that was generated from the dependent job 
# submission script and goes through and cancels all of the jobs.
#
# Note: This script was written to run on LLNL's Vulcan computer. It may work on other
# machines with similar job scheduler commands.
#
#Author: Christopher Neal
#
# Date: 02-18-2015
# Updated: 03-26-2015
#
#######################################################################################

# This script should be executed in the directory that contains the jobID.txt file that
# was generated from the dependent job submission script.

echo 'Running Dependent Job Cancellation Script'


#Loop through the file and read each line one-at-a-time and cancel the corresponding job.
while read line
do

    canceljob $line   

done < jobIDs.txt 

#Remove the file containing the list of job submission IDs
rm -rf jobIDs.txt

echo 'Program Finished...Ending'


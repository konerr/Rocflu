#!/bin/bash
#######################################################################################
# The purpose of this script is to submit several depdendent jobs. The job submission
# script will be the same; this just resubmits it several time for code restarting
# purposes.
#
#
#Author: Christopher Neal
#
# Date: 02-11-2015
# Updated: 02-18-2015
#
#######################################################################################

#This script should be executed in the directory that contains the job submission script.
# It takes 2 arguments:
# 1.)The name of the job submission script that you want to submit
# 2.) The number of resubmissions that you want i.e. an integer 
#
# The code creates a file that contains all of the job IDs from the jobs that are submitted.
#
echo 'Running Dependent Job Submission Script'

echo "Submitting Job File: $1"

for (( i=1; i<=$2; i++))
do

  echo "Submitting Job: $i"

  if [ "$i" -eq "1" ]; then
    
    msub $1 | tee submission_temp.txt #Copy output from terminal into file(COPIES)
    one=$(awk '/./{line=$0} END{print line}' submission_temp.txt) #Extract last numeric line of temp file
    #echo $one #For debugging
    
    echo $one > jobIDs.txt #Write job ID to file

  else
    
    msub -l depend=afterany:$one $1 | tee submission_temp.txt
    two=$(awk '/./{line=$0} END{print line}' submission_temp.txt)
    #echo $two #For debugging

    one=$two

    echo $one >> jobIDs.txt #Append job ID to file

  fi

done


#Clean up temporary files
rm -rf submission_temp.txt

echo 'Program Finished...Ending'


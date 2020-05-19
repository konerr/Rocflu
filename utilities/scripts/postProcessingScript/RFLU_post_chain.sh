#!/bin/bash
#######################################################################################
# Purpose: The purpose of this script is to submit several independent post-processing jobs. 
#          The job submission script will be the slighly the same. The only difference
#          will be the timestamp that the executable will be processing.
#
#
# Author: Christopher Neal
#
# Date: 03-27-2015
# Updated: 03-27-2015
#
#######################################################################################

# This script should be executed in the directory that contains the template job 
# submission script and a file called Extracted_Data_Labels.txt which contains the
# the timesteps that are to be processed.
#
# It takes the following 2 arguments:
# 1.) The solution stride that you want to post-process
# 2.) The integer value of the solution that you want to start processing from
#	i.e. if argument 1 is 1 and argument 2 is 6 then the script will process
#	all solutions INCLUDING after the 6th one in the directory. This feature is 
#	useful for processing additional solutions without re-processing every solution 
#	file in the working directory.
#
echo 'Running Rocflu Post-Processing Chain Submission Script'


##Extract Solution Times####
echo "Extracting List of Solution Times From .mixt Solution Files"

#Put all files in current directory into a text file
ls *mixt.cva*00001* >temp.txt


#Trim filename up the first underscore.
sed -n -i "s/[^_]*_//p" temp.txt


#Trim filename up to the second underscore.
sed -n -i "s/[^_]*_//p" temp.txt


#At this point we have all of the times, but they are not sorted.
#Sort the times that are in Extracted_Data_Labels.txt

sort -k1g temp.txt > Extracted_Data_Labels.txt

#Remove temporary files
rm -rf temp.txt

###


# Copy the job_RFLU_postchain.msub to a temporary file to prevent the template
# from being permanently changed.
cp -rf  job_RFLU_postchain.msub tempScript.msub

#Loop over times in the file Extracted_Data_Labels.txt and submit jobs
c=1
while read line1
do

if [ $(($c%$1)) == 0 ]; then #Check to see if the file is a multiple of the stride # input

   if [ -z "$2" ]; then #If argument returns true is the second input argument is empty

      echo "Submitting job file $c for post-processing time: $line1 "
   
      #Replace job log file line in submission script with unique identifier
      ReplaceString="#MSUB -o job_RFLU_postchain$c.log"
      sed -i "s/#MSUB -o.*/$ReplaceString/" job_RFLU_postchain.msub


      # The variable $line1 contains the timestamp that we want to run rflupost for
      # Replace the TimeStamp entry in the job submission template with the actual time
      # and submit job

      ReplaceString="$line1"  

      sed -i "s/TimeStamp.*/$ReplaceString/" job_RFLU_postchain.msub

      #Submit the job file
      msub job_RFLU_postchain.msub

      #Copy over the altered job file with the temporary template tempScript.msub
      cp -rf  tempScript.msub job_RFLU_postchain.msub
  
   else
   
      if [ $c -ge $2 ]; then  #Process files starting from file number stored in $2
         echo "Submitting job file $c for post-processing time: $line1 "

         #Replace job log file line in submission script with unique identifier
         ReplaceString="#MSUB -o job_RFLU_postchain$c.log"
         sed -i "s/#MSUB -o.*/$ReplaceString/" job_RFLU_postchain.msub


         # The variable $line1 contains the timestamp that we want to run rflupost for
         # Replace the TimeStamp entry in the job submission template with the actual time
         # and submit job

         ReplaceString="$line1"

         sed -i "s/TimeStamp.*/$ReplaceString/" job_RFLU_postchain.msub

         #Submit the job file
         msub job_RFLU_postchain.msub

         #Copy over the altered job file with the temporary template tempScript.msub
         cp -rf  tempScript.msub job_RFLU_postchain.msub
      fi #$c -ge $2

   fi #-z "$2"

fi #$(($c%$1)) == 0 

#Increment the counting variable, c
c=$(($c+1))

done < Extracted_Data_Labels.txt

# Write over the altered job_RFLU_postchain.msub file with the unaltered version
# stored in tempScript.msub

cp -rf  tempScript.msub job_RFLU_postchain.msub 

# Delete temporary file
rm -rf tempScript.msub

echo 'Program Finished...Ending'


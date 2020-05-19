#!/bin/bash
#######################################################################################
# The purpose of this script is to delete solution files in a directory. The script
# looks at the solutions in the directory and performs a strided deletion of the files.
# For example: If the argument to the script is 2 then it will delete every second file.
#
#
# Author: Christopher Neal
#
# Date: 06-17-2015
# Updated: 06-23-2015
#
#######################################################################################

#This script should be executed in the directory that contains the solution files.
# It takes 2 arguments:
# 1.) The stride of the deletion in the form of an integer 
# 2.) The casename of the solution i.e. cylds.rin  pass the function cylds
#
# Example Use: RFLU_delete.sh 2 would do the following to the solutions in a directory
#
#	Original Solutions			New Solutions
#	0.00000E+00				0.00000E+00
#	1.00000E-01				2.00000E-01
#	2.00000E-01				4.00000E-01
#	3.00000E-01				6.00000E-01
#	4.00000E-01
#	5.00000E-01
#	6.00000E-01
#
# If every solution time is sorted from low to high and assigned a rank that starts at 1 and
# increases by 1, then it is easy to see that the code simply deletes times that have a rank that 
# is divisible by the argument that is passed to the script.

echo 'Running Rocflu Solution Deletion Script'

echo "Removing a solution for every $1 solutions"

#Generate the extracted data labels
echo "Extracting Times From .mixt Solution Files"

#Put all files in current directory into a text file
for f in *; do echo "$f"; done >temp.txt

#Print the occurences of the first processor solution files to a new file
sed -n "/.mixt.cva_00001_/p" temp.txt >temp2.txt


#Trim filename up the first underscore.
sed -n -i "s/[^_]*_//p" temp2.txt

#Trim filename up to the second underscore.
sed -n -i "s/[^_]*_//p" temp2.txt

#Remove temporary files
rm -rf temp.txt

#At this point we have all of the times, but they are not sorted.
#Sort the times that are in Extracted_Data_Labels.txt
sort -k1g temp2.txt > Extracted_Data_Labels.txt

#Remove temporary files
rm -rf temp2.txt


#Proceed to stage where code deletes solution files
c=1
while read line1
do

  #Debug
  #echo "Variable is: $(($c%$1)) :  Count is $c"

  if ! [ -z "$line1" ] && [ $(($c%$1)) == 0 ]; then #check if file is multiple of the stride input number

	#If there is a match we should not copy this into the new file 
	echo "Deleting line $c From .rin file & Solutions with Time: $line1"

	#Delete all files with the matching timestamp
        rm -rf *$line1*


  else

	if [ $c == 1 ]; then

		#Copy this line into a new file
	        sed  "${c}q;d" $2.rin > temp_times.txt

	else

		#Copy this line into a new file
        	sed  "${c}q;d" $2.rin >> temp_times.txt

	fi

  fi

#Increment counter
c=$(($c+1))

done<Extracted_Data_Labels.txt

#Copy over the temporary file to the original .rin file
mv temp_times.txt $2.rin

#Clean up temporary files
rm -rf temp_times.txt

echo 'Program Finished...Ending'


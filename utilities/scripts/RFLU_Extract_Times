#!/bin/bash

################################################################################
#
# Purpose: To extract the time stamps on all of the .mixt.cvs solution files in a
#	   directory and store them sorted from smallest to largest into a file
#          called Extracted_Data_Labels.txt.
#
#
#
#
# Author: Christopher Neal
# Date:	  1/15/2015
# Updated: 3/12/2015
################################################################################

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


echo "Script Finished Running. Check Extracted_Data_Labels.txt for list of timestamps."
#-----------------------------------------------------------------------------

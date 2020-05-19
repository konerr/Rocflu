#! /usr/bin/env python

# Purpose:	Reads in and change the names of the rocflu post Ensight Files for use of 
#		Paraview.
#
# Description:  RocFlu post can output Ensight files, however the case file and filenames of 
#		the post output files are in an incorrect format for Paraview to open all the 
#		output files at once. The program rectifies this by doing two important steps. 
#		
#		Step 1: Search through the files, of the current working directory, for the 
#			post output files and change the names of the files to a proper format
#			cylds."property"_00001_"timestep" => cylds."property"_"intstep" 
#			for an explination of intstep see Notes(3)
#		
#		Step 2: Build a new case file with the proper syntax so Paraview can open all
#			the timesteps from just the case file. 

# Input: 	None.

# Output: 	None.

# Notes:(1)	Since this program is designed to rename the RocFlu Ensight Gold output format
#		it is thus assumed that the script, which is called "updateRFLUEnsightOutput", 
#		was executed to generate the output files. If the script was not called the  
#		the program will not work. The bash script also outputs a text file called 
#		Extracted_Data_Labels.txt. This file is crucial for the renaming of the files.
#		
#	(2)	This program and the bash script assumes that RFLU Post will output the post
#		files in the form of Ensight Gold format. In cylds.inp the #POST flag section
#		there must be an OUTFLOW flag with a value of "2", else there will be no files
#		to rename.
#
#	(3)	The file format that RocFlu post outputs is in the form of 
#		cylds."property"_00001_"timestep", where property is a scalar quantity, like 
#		Desnity or Pressure, and timestep is the timestep of that solution file. The 
#		timestep in the name has been truncated by the bash script. 
#		
#	(4)	Step 1 assumes that the file output of the RocFlu post is in the form of 
#               cylds."scalar"_00001_"timestep" if the output filename is not in that form the
#               program will not find the post output files
#		
#	(5)	Paraview will not recognize these filenames as a set of timesteps and will only
#		be able to open one at a time. Paraview will recoginze files as a series of
#		timesteps if there is an accending order of integers in the solution filenames
#		with "buffer" zeros in front of the integer. Example: 0001, 0002, 0003,...,0010
#		0011,...,0100,0101,0102. The buffer zeros should be in place for at least the 
#		maximum	amount of timesteps that will be needed. Example: If 900 timesteps,
#		the first time step should be 001, 0001,  or 00001, but not 01. 
#	
#	(6)	Step 2 Assumes that the data in Extracted_data_Labels.txt is sorted from 
#		smallest to largest. At the time of this program writing this should be the
#		case. This sort assumption is crucial for determing in which order Paraview 
#		should open up each solution file. The index location of the timestep in the
#		time step array is the order in which Paraview will open the files. Starting 
#		at zero. Example: time_step[0] is the first timestep, time_step[1] is the
#		second time step. 
#
#	(7)	Step 2 rebuilds the case file by overwriting over the old one. In building the
#		case file it is assumed that the geometry of the run does not change and that  
#		there is only one time set.
#
#	Authors: Christopher Neal, Bradford Durant
#	Date   : 10/1/2014
#	Updated: 1/15/2015
#
########################################################################

####### ASSUMPTIONS #######

####### updateRFLUEnsightOutput was excecuted 

####### Assumes that rflupost outputs files in the form of <casename>."scalar"_00001_"timestep"
####### Example: cylds.r_00001_0.00000E+00

####### Extracted_Data_Labels.txt exists and is sorted from smallest to largest

####### In building the case file:      
####### Assumes that geometry of the run does not change
####### Assumes that there is only one time sets

#for use of reading and writing files
import os
 
#Read in the file with the time step times 
timeSteps_filename = "Extracted_Data_Labels.txt"

#Declaring varables
time_step = []
step_zeros = '000'
Casename = "ctint"

#If file exists, open it
if timeSteps_filename in os.listdir(os.getcwd()):
	input_file = open(timeSteps_filename)
	
	#Go through file and put time steps into an array
	#Assumed the times are already a sorted from bash script
	#Grabbing only the time value not the "\n" value	
	for line in input_file:
		
		#Grabbing only the time value not the "\n" value 
		time_step.append(line[:11])
	input_file.close()
	#print (time_step)
else:
	print("Extracted_Data_Labels.txt not found.")
	print("bash script updateRFLUEnsightOutput did not run as assumed.") 

#When RFLU post creates Ensight Gold files it creates a case file.
#If case file is not found then RFLU did not output as assumed
if os.path.exists("./"+Casename+".case_00001") == False:
	print("Original output case file not found")
        print("Either updateRFLUEnsightOutput was not executed as assumed")
	print("or RFLU may have not output in EnSight Gold format")
        print("Check flag in #Post in <casename>.inp")
        print("Flag should be OUTFORMAT 2")

Old_FileName = Casename+".case_00001"
New_FileName = Casename+".case.00001"

os.rename(Old_FileName,New_FileName)

####### Assumes that rflupost outputs files in the form of <casename>."scalar"_00001_"timestep"
####### Example: cylds.r_00001_0.00000E+00

####### Coverts file name to <casename>."property"_"intstep"
Schlieren_data = False # Sometimes extra data is output for 2D simulations

base_len = len(Casename+".r_00001_") # smallest string length computed using casename
case_len = len(Casename) # Length of casename

#Go through files in current directory
for filename in os.listdir(os.getcwd()):
	
	#Set flag to default
	match_found = False

	#If the scalar data files has default name from RFLUPost
	if filename.startswith(Casename+".r_00001_"):

		#Grab the tail end of the file name, which is the time step time
		time_check = filename[base_len:]
		match_found = True
		
		#Grab the front part of the file name up till _00001
		front_filename = filename[:case_len+2]

	elif filename.startswith(Casename+".rv_00001_"):
		time_check = filename[base_len+1:]
		match_found = True
		front_filename = filename[:case_len+3]

	elif filename.startswith(Casename+".rE_00001_"):
	 	time_check = filename[base_len+1:]
                match_found = True
                front_filename = filename[:case_len+3]

	elif filename.startswith(Casename+".p_00001_"):
		time_check = filename[base_len:]
		match_found = True
		front_filename = filename[:case_len+2]
	
	elif filename.startswith(Casename+".T_00001_"):
		time_check = filename[base_len:]
		match_found = True
		front_filename = filename[:case_len+2]
	
	elif filename.startswith(Casename+".a_00001_"):
		time_check = filename[base_len:]
		match_found = True
		front_filename = filename[:case_len+2]
	elif filename.startswith(Casename+".schl_00001_"):
		time_check = filename[base_len+3:]
		match_found = True
		front_filename = filename[:case_len+5]
		Schlieren_data = True
	elif filename.startswith(Casename+".shad_00001_"):
                time_check = filename[base_len+3:]
                match_found = True
                front_filename = filename[:case_len+5]
                Schlieren_data = True
	elif filename.startswith(Casename+".intf_00001_"):
                time_check = filename[base_len+3:]
                match_found = True
                front_filename = filename[:case_len+5]
                Schlieren_data = True
	else:
		match_found = False

	#If we found a file with the inital rocflu post name set up the buffer zeros for approperiate counting 
	if match_found == True:
		count = 0
		time_matchFound = False
		
		#Loop through the time_step array
		while count <= (len(time_step)-1) and time_matchFound == False:
			
			#If the timestep in the filename matches the timestep in Extracted_Data_Labels.txt
			#Apply the approperiate buffer zeros
			if time_step[count] == time_check:
				if count > 999:
					mid_filename = ''
				elif count > 99:
					mid_filename = step_zeros[:1]
				elif count > 9:
					mid_filename = step_zeros[:2]
				else:
					mid_filename = step_zeros

				time_matchFound = True
				
				#The position the time step was found in the array from Extracted_Data_Labels.txt
				#determines which iteration of the solution it is. 
				count_s = str(count)
				
				#Build the new filename in the form of cylds."property"_"intstep"
				full_filename = front_filename + '_'  + mid_filename  + count_s # + '_' + time_check
				os.rename(filename,full_filename)

			count +=1
	
	#If file name found but its time step did not match from Extracted_Data_Label.txt
	if match_found == True and time_matchFound == False:
		print("Error: time step not found!")
		print("There may be a solution file that has not been renamed")
		print("And thus may not be opened by Paraview")

#Rebuild the case file so that paraview can open as a multi file transient set
#The name of the case file that Paraview will recognize is seperated by "."
case_filename = Casename+".case.00001"

#The file format is hard coded in this program.

####### In building the case file:	
####### Assumes that geometry of the run does not change
####### Assumes that there is only one time set

lines_2add = []

if Schlieren_data == False:
	lines_2add.append('FORMAT')
	lines_2add.append('type: ensight gold')
	lines_2add.append('GEOMETRY')
	lines_2add.append('model: ./'+Casename+'.geo_00001_000000')
	lines_2add.append('VARIABLE')
	lines_2add.append('scalar per element: Density ./'+Casename+'.r_****')
	lines_2add.append('scalar per element: Energy ./'+Casename+'.rE_****')
	lines_2add.append('scalar per element: Pressure ./'+Casename+'.p_****')
	lines_2add.append('scalar per element: Temperature ./'+Casename+'.T_****')
	lines_2add.append('scalar per element: Soundspeed ./'+Casename+'.a_****')
	lines_2add.append('vector per element: Momentum ./'+Casename+'.rv_****')
	lines_2add.append('TIME')
	lines_2add.append('time set:             1')
	lines_2add.append('number of steps:      ' + str(len(time_step)))
	lines_2add.append('filename start number:  0')
	lines_2add.append('filename increment:   1')
	lines_2add.append('time values:  ') 
else:
	lines_2add.append('FORMAT')
        lines_2add.append('type: ensight gold')
        lines_2add.append('GEOMETRY')
        lines_2add.append('model: ./'+Casename+'.geo_00001_000000')
        lines_2add.append('VARIABLE')
	lines_2add.append('scalar per element: Density ./'+Casename+'.r_****')
        lines_2add.append('scalar per element: Energy ./'+Casename+'.rE_****')
        lines_2add.append('scalar per element: Pressure ./'+Casename+'.p_****')
        lines_2add.append('scalar per element: Temperature ./'+Casename+'.T_****')
        lines_2add.append('scalar per element: Soundspeed ./'+Casename+'.a_****')
        lines_2add.append('vector per element: Momentum ./'+Casename+'.rv_****')
	lines_2add.append('scalar per element: Interferometry ./'+Casename+'.intf_****')
	lines_2add.append('scalar per element: ShadowGraph ./'+Casename+'.shad_****')
	lines_2add.append('scalar per element: Schlieren ./'+Casename+'.schl_****')
        lines_2add.append('TIME')
        lines_2add.append('time set:             1')
        lines_2add.append('number of steps:      ' + str(len(time_step)))
        lines_2add.append('filename start number:  0')
        lines_2add.append('filename increment:   1')
        lines_2add.append('time values:  ')


#After "time values:" the time steps need a numerical name to associate which timestep is what time.
# Example time step 5.123E-02 can be called time step 0.05123  

#Open the file, 'w+' means to overwrite any exisiting files by that name
write_file = open(case_filename,'w+')
write_count = 0

#Write the file as hard coded above.
while write_count <= (len(lines_2add)-1):
	write_file.write(lines_2add[write_count])
	write_count +=1
	write_file.write("\n")

#Add the timesteps' numerical names after the line 'time values:' in the case file 
count = 0
while count <= (len(time_step)-1):
	write_file.write(time_step[count])
	write_file.write("\n")
	count +=1

write_file.close()

print("\n Program has finished... \n")


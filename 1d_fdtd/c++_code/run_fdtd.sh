#!/usr/bin/env bash

#This will only work in the release build!!!
#Sets the directory of the output file as well as the output file
fdtd_out=main.out
fdtd_bin=./src/bin/$fdtd_out
echo ================================================================
echo Directory of compiled source code: $fdtd_bin

#Set the directory and filename of the input file
input_file=inputs.csv
input_dir=../inputs/$input_file
#input_dir=../../../inputs/$input_file
#input_dir=manual

echo Directory of input file: $input_dir

#Set the output directory
output_dir=../../../outputs/

#The code here has no exception handling so be careful!!
echo ================================================================
#Boundary condition: either 'pabc' or 'dirichlet'
read -p 'Select the Boundary Condition [pabc/dirichlet]: ' boundary_cond

#Source Excitation method: 'hard', 'soft', or 'tfsf'
read -p 'Select Source Excitation Method [hard/soft/tfsf]: ' excitation_method

#Custom name for output file: any string of char
read -p 'Enter custom name for output file/s: ' custom_name

#Output file type: 'csv' or  'npy'
read -p 'Select output file type [csv/npy/hdf5]: ' output_file_type

#Type of algorithm to be used: 'fdtd' or 'fdtd-schwarz'
read -p 'Select algorithm [fdtd/fdtd-schwarz]: ' algorithm

#Number of subdomains for 'fdtd-schwarz' algorithm [greater than or equal to 1]
read -p 'How many subdomains will be used: ' num_subdomains

#Toggle if multithreading will be used 
read -p 'Will multithreading be used? [true/false]: ' multithreading

#Enter overlap size (must be greater than 0 and less than Nz)
read -p 'How much overlap between subdomains? [greater than 0]: ' overlap_size

#Check if it will be comprehensive in saving data or not
read -p 'Comprehensive data saving? [true/false]: ' comprehensive

#Username or Name of user
read -p 'Username/Name: ' user_name

#Short description of the simulation
read -p 'Description of the simulation: ' sim_description

echo ================================================================
echo Running the compiled program.....


#Execute the program
#Format of input arguments: (1) input file directory (2) Boundary Condition (3) Excitation method (4) Custom Output Filename (5) Output File type (6) Output directory
log_file=$(date '+%Y-%m-%d')_$custom_name.log

$fdtd_bin $input_dir $boundary_cond $excitation_method $custom_name $output_file_type $output_dir $algorithm $num_subdomains $multithreading $overlap_size $comprehensive $user_name $sim_description |& tee $output_dir/$log_file
echo Check the simulation logs in $output_dir/$log_file to check the simulation details
echo ================================================================

#Execute the python script
#Get the input arguments; Format: (1) input directory (2) Custom file name (3) Output file type (4)Source Name (for image) (5) Output directory (for plots)
python_file=plotly_subplot.py
py_output_dir=$output_dir/plots/
#echo ================================================================
#read -p "Enter a custom name for the source image: " source_name
echo ================================================================
echo Running plotting script....
py_input=$(date '+%Y-%m-%d')_$custom_name.hdf5
python3 $python_file $output_dir $py_input #$output_file_type $py_output_dir |& tee -a $output_dir/$log_file
echo Check the simulation logs in $output_dir/$log_file to check the plotting details
#echo ================================================================
#read -p "Enter output filename for video file: " video_name
##echo Running ffmpeg...
#ffmpeg -f image2 -framerate 200 -i $output_dir/plots/E_H_FFT_images_%07d.png -s 1920x1080 $output_dir/$video_name.mp4
#echo ================================================================
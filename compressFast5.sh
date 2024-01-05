#!/bin/bash


fast5=$1
alphaLocation=$2

# for reads in "$fast5"/*/*/*.fast5;

# 	do 
# 		cp "$reads" "$fast5"
		
# 	done 




# Set the directory path to the folder containing the files to be compressed
# DIR_PATH=$fast5
 
# # Set the name prefix of the output archive files
# ARCHIVE_PREFIX="archive"
 
# # Set the maximum number of files per batch
# MAX_FILES=50
 
# # Change directory to the specified path
# cd "$DIR_PATH"
 
# # Get a list of all files in the directory
# files=( * )
 
# # Calculate the number of batches of files
# num_batches=$(( (${#files[@]} + $MAX_FILES - 1) / $MAX_FILES ))
 
# # Loop through each batch of files
# for (( i=0; i<$num_batches; i++ )); do
#     # Set the start and end indices of the batch
#     start=$(( $i * $MAX_FILES ))
#     end=$(( ($i + 1) * $MAX_FILES - 1 ))
     
#     # Check if the end index exceeds the number of files
#     if (( $end >= ${#files[@]} )); then
#         end=$(( ${#files[@]} - 1 ))
#     fi
     
#     # Create a compressed archive file for the batch of files
#     archive_name="${ARCHIVE_PREFIX}_${i}.zip"
#     zip "$archive_name" "${files[@]:$start:$MAX_FILES}"
# done


for readArxv in "$fast5"*.zip;

	do

		scp "$readArxv" dd87a@alpha2.cvr.gla.ac.uk:/home3/dd87a/"$alphaLocation"

	done







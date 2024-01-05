#!bin/bash 

for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ; 

do

	folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

	for file in "$folder"*.fastq ; 

		do
			NanoStat --fastq *.fastq --outdir statReports -n "$folderName".txt

		done 

done
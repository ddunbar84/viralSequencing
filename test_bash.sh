#!bin/bash 

for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ; 

do

echo "$folder"

	for file in "$folder"*.fastq ; 

	do

 	echo "$file"

	done 

done



#for d in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/*.fastq ; do

	#for d in /home3/dd87a/{fast5}/demulti_trim/*/*.fastq ; 
		#do

			 #NanoStat --fastq *.fastq --outdir statReports -n barcode01.txt

			#echo "$d" # test that for loop works

		#done 


for file in /home3/dd87a/FHV_iso_26may23/*Chop* ;

do

        fileName=`echo "$file" | sed -e 's/^.*\(23\)\///'| sed -e 's/\.fastq$//'`

        NanoStat --fastq "$file" --outdir /home3/dd87a/FHV_iso_26may23/statreport -n "$fileName".txt

done
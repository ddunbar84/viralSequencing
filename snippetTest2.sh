
# Stout to file for error checking

{
##### tmux
##### run process
##### ctrl +b, d (detach from tmux)

##### in shell tmux attach to get back into tmux screen

##### Log in to Alpha2
##### screen -S GPU3
##### ssh gpu3 (prompts for host then password)
##### nvidia-smi

project=$1 # project should be the project folder (and all files should be placed in the subfolder "allFast5")
refseq=$2 # Refseqs for project should be stored in the project folder level
# $3 = input for barcoder
# $4 = output for barcoder (input for porechop?)
# $5 = quality trimmed reads input for minimap?
# $6 = nanopolish??
# $7 =
# $8 = 
# $9 =  

##### ACTIVATE "ont" CONDA ENV - NEEDED FOR DATA PROCESSING STEPS
source /home3/dd87a/miniconda3/bin/activate ont

mkdir /home3/dd87a/"$project"/basecalled;
mkdir /home3/dd87a/"$project"/demulti_trim;
mkdir /home3/dd87a/"$project"/assemblies;
mkdir /home3/dd87a/"$project"/draftAssembly;
mkdir /home3/dd87a/"$project"/medaka/annotatedVCF;


################ PASTE SNIPPET HERE ################

for filename in /home3/dd87a/"$project"/medaka/*VC/????_medaka.annotated.vcf;

do
	
	vcfName=`echo "$filename" | sed -e 's/^.*VC\///'`

	cp "$filename" /home3/dd87a/"$project"/medaka/annotatedVCF/"$vcfName"

done


} > snippetTestLog.txt




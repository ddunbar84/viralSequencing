#tmux
#run process
#ctrl +b, d (detach from tmux)

#in shell tmux attach to get back into tmux screen

#Log in to Alpha2
#screen -S GPU3
#ssh gpu3 (prompts for host then password)
#nvidia-smi

project=$1 # project should be the project folder (and all files should be placed in the subfolder "allFast5")
refseq=$2 # Refseqs for project should be stored in the project folder level
# $3 = input for barcoder
# $4 = output for barcoder (input for porechop?)
# $5 = quality trimmed reads input for minimap?
# $6 = nanopolish??
# $7 =
# $8 = 
# $9 =  

### ACTIVATE "ont" CONDA ENV - NEEDED FOR DATA PROCESSING STEPS
source /home3/dd87a/miniconda3/bin/activate ont

mkdir /home3/dd87a/"$project"/basecalled;
mkdir /home3/dd87a/"$project"/demulti_trim;
mkdir /home3/dd87a/"$project"/assemblies;

### BASECALLING ###

# BASECALL RNA ON CPU
# guppy_basecaller --input_path /home3/dd87a/"$project"/ --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109

# # BASECALL RNA ON GPU
# /software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/"$project"/ --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109 -x "cuda:0"

# BASECALL RNA ON CPU
# guppy_basecaller --input_path /home3/dd87a/"$project"/ --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96

# BASECALL ON GPU
# /software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/"$project"/allFast5 --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96 -x "cuda:0"

#cat file*.fastq > bigfile.fastq

# Concatenate all fastq basecalled files into one file to feed into nanopolish
cat /home3/dd87a/"$project"/basecalled/pass/fastq*.fastq > basecalled.fastq

# Nanopolish raw data to try and improve read quality 

nanopolish index -d /home3/dd87a/"$project"/allFast5 -f /home3/dd87a/"$project"/basecalled/ /home3/dd87a/"$project"/basecalled/pass/basecalled.fastq
#nanopolish index -d /home3/dd87a/FHV_scriptTest/allFast5 -f /home3/dd87a/FHV_scriptTest/basecalled/ /home3/dd87a/FHV_scriptTest/basecalled/pass/basecalled.fastq

# ### DEMULTIPLEXING ###

# # DEMULTIPLEX BASECALLED READS
guppy_barcoder -r --input_path /home3/dd87a/"$project"/basecalled/pass --save_path /home3/dd87a/"$project"/demulti_trim --trim_barcodes 
# #guppy_barcoder -r --input_path /home3/dd87a/FHV_test_28feb23/basecalled/pass --save_path /home3/dd87a/FHV_test_28feb23/basecalled/pass/demultiplexed_trimmed --trim_barcodes 

# ### INITIAL QC STEP FOR EACH BARCODE ###

# RUN QC REPORT FOR EACH BARCODE
 for folder in /home3/dd87a/"$project"/demulti_trim/*/ 
# #for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

 do

         folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

         NanoStat --fastq "$folder"*.fastq --outdir /home3/dd87a/"$project"/statreport -n "$folderName".txt

 done

### ASSESS READ QUALITY AND TRIM POOR QUALITY READS

## change this to make a folder for the chop data

# Use porechop to assess read quality and trim the poor quality reads
# requires input dir, output file and flag to discard reads with adapters in the middle (this is required for nanopolish)

#for folder in /home3/dd87a/"$project"/demulti_trim/*/ 
for folder in /home3/dd87a/"$project"/demulti_trim/*/ ;

do

       folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

       /software/Porechop-v0.2.4/porechop-runner.py -i "$folder" -o /home3/dd87a/"$project"/demulti_trim/"$folderName"/"$folderName"Chop.fastq --discard_middle

done


### RUN QC REPORT FOR EACH BARCODE POST PORECHOP

for sample in /home3/dd87a/"$project"/demulti_trim/*/*Chop* ; 
#for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

do

        folderName=`echo "$sample" | sed -e 's/^.*\(barcode*\/\)\///'| sed -e 's/\.fastq$//'`

        NanoStat --fastq "$sample" --outdir /home3/dd87a/"$project"/statreport -n "$folderName".txt

done

# POSSIBLY ADD A CHOPPER SECTION IN HERE TO FILTER THE READS TO LONGER READS (>300bp) ???

# FLYE ASSEMBLY TO USE WITH MEDAKA


# CANU ASSEMBLY MAY PROVIDE DIFFERENT RESULTS WITH MEDAKA (also canu is not currently working??)
# canu -d /home3/dd87a/FHV_scriptTest/canu -p canu-assembly useGrid=false -nanopore-raw /home3/dd87a/FHV_scriptTest/basecalled/pass/basecalled.fastq  genomeSize=135k minReadLength=300 minOverlapLength=20 stopOnReadQuality=false useGrid=false

#need to add full ref to 
quast.py <input files> -t 8 -o <outdir> -r <refgenome> <files_with_contigs>

# use minimap2 (or flye) to map reads to reference genomev1rology

# -ax CIGAR and output to SAM, (x) map-ont - map ONT data, align noisy long reads 10% error rate to ref
# -bS output to BAM, compatibitity check          -F exclude FLAG'd alignments mapping reads
# -o output file 

for folder in /home3/dd87a/"$project"/demulti_trim/*/;

do

       folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

       minimap2 -ax map-ont "$project"/"$refseq" "$folderName"Chop.fastq | samtools view -bS -F 4 | samtools sort -o /home3/dd87a/"$project"/assemblies/"$folderName"draft.bam

done

#samtools view *.bam | head

### RUN QC REPORT FOR EACH BARCODE POST ASSEMBLY

for algn in /home3/dd87a/"$project"/assemblies/*.bam ; 
#for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

do

        algnName=`echo "$algn" | sed -e 's/\.bam$//'`
        
        NanoStat --bam "$algn" --outdir /home3/dd87a/"$project"/statreport -n "$algnName".txt

done


# # Add a WIMP mapping protocol?

# # Add in a QC check step?

# #nanopolish-v0.14.0  

# #/software/nanopolish-v0.14.0/nanopolish


#nanopolish vcf2fasta -g draft.fa polished.*.vcf > polished_genome.fa

# ivar call consensus / variant calling?


#!bin/bash 

#tmux
#run process
#ctrl +b, d (detach from tmux)

#in shell tmux attach to get back into tmux screen

#Log in to Alpha2
#screen -S GPU3
#ssh gpu3 (prompts for host then password)
#nvidia-smi


# $1 = input path for basecaller
fast5 = $1
refseq = $2
# $3 = input for barcoder
# $4 = output for barcoder (input for porechop?)
# $5 = quality trimmed reads input for minimap?
# $6 = nanopolish??
# $7 =
# $8 = 
# $9 =  

### ACTIVATE "ont" CONDA ENV - NEEDED FOR DATA PROCESSING STEPS
source /home3/dd87a/miniconda3/bin/activate ont


### BASECALLING ###

# BASECALL RNA ON CPU
# guppy_basecaller --input_path /home3/dd87a/{fast5}/ --save_path /home3/dd87a/{fast5}/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109

# # BASECALL RNA ON GPU
# /software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/{fast5}/ --save_path /home3/dd87a/{fast5}/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109 -x "cuda:0"

# BASECALL RNA ON CPU
# guppy_basecaller --input_path /home3/dd87a/{fast5}/ --save_path /home3/dd87a/{fast5}/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96

# BASECALL ON GPU
/software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/{fast5}/ --save_path /home3/dd87a/{fast5}/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96 -x "cuda:0"


### DEMULTIPLEXING ###

# DEMULTIPLEX BASECALLED READS
guppy_barcoder -r --input_path /home3/dd87a/{fast5}/basecalled/pass --save_path /home3/dd87a/{fast5}/demulti_trim --trim_barcodes 
#guppy_barcoder -r --input_path /home3/dd87a/FHV_test_28feb23/basecalled/pass --save_path /home3/dd87a/FHV_test_28feb23/basecalled/pass/demultiplexed_trimmed --trim_barcodes 

### INITIAL QC STEP FOR EACH BARCODE ###

# RUN QC REPORT FOR EACH BARCODE
for folder in /home3/dd87a/{fast5}/demulti_trim/*/ 
#for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

do

        folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

        NanoStat --fastq "$folder"*.fastq --outdir /home3/dd87a/FHV_iso_26may23/statreport -n "$folderName".txt

done

### ASSESS READ QUALITY AND TRIM POOR QUALITY READS

# Use porechop to assess read quality and trim the poor quality reads
# requires input dir, output file and flag to discard reads with adapters in the middle (this is required for nanopolish)

#for folder in /home3/dd87a/{fast5}/demulti_trim/*/ 
for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

do

        folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

        /software/Porechop-v0.2.4/porechop-runner.py -i "$folder" -o /home3/dd87a/{fast5}/demulti_trim/"$folderName"/"$folderName"Chop.fastq --discard_middle

done

### RUN QC REPORT FOR EACH BARCODE POST PORECHOP

for folder in /home3/dd87a/{fast5}/demulti_trim/*/ 
#for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

do

        folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

        NanoStat --fastq "$folder"Chop.fastq --outdir /home3/dd87a/FHV_iso_26may23/statreport -n "$folderName"Chop.txt

done

# use minimap2 (or flye) to map reads to reference genome
# -ax CIGAR and output to SAM, (x) map-ont - map ONT data, align noisy long reads 10% error rate to ref
# -bS output to BAM, compatibitity check          -F exclude FLAG'd alignments mapping reads
# -o output file 

for folder in /home3/dd87a/{fast5}/demulti_trim/*/ 
#for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

do

        folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

        minimap2 -ax map-ont {refseq} "$folderName"Chop.fastq | samtools view -bS -F 4 | samtools sort -o draftGenome.bam

done
# Add a WIMP mapping protocol?

# Add in a QC check step?

nanopolish-v0.14.0  

/software/nanopolish-v0.14.0/nanopolish

nanopolish vcf2fasta -g draft.fa polished.*.vcf > polished_genome.fa

# ivar call consensus / variant calling?




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


##### ### BASECALLING ###

##### BASECALL RNA ON CPU
# guppy_basecaller --input_path /home3/dd87a/"$project"/ --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109

##### BASECALL RNA ON GPU
# /software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/"$project"/ --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109 -x "cuda:0"

##### BASECALL RNA ON CPU
# guppy_basecaller --input_path /home3/dd87a/"$project"/ --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96

##### BASECALL ON GPU
# /software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/"$project"/allFast5 --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96 -x "cuda:0"


# ##### ### DEMULTIPLEXING ###

# ##### DEMULTIPLEX BASECALLED READS

# guppy_barcoder -r --input_path /home3/dd87a/"$project"/basecalled/pass --save_path /home3/dd87a/"$project"/demulti_trim --trim_barcodes 
# ##### # guppy_barcoder -r --input_path /home3/dd87a/FHV_test_28feb23/basecalled/pass --save_path /home3/dd87a/FHV_test_28feb23/basecalled/pass/demultiplexed_trimmed --trim_barcodes 

# ##### ### INITIAL QC STEP FOR EACH BARCODE ###

# ##### RUN QC REPORT FOR EACH BARCODE

#  for folder in /home3/dd87a/"$project"/demulti_trim/*/ 
# ### for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

#  do

#          folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

#          NanoStat --fastq "$folder"*.fastq --outdir /home3/dd87a/"$project"/statreport -n "$folderName".txt

#  done

# ##### ### ASSESS READ QUALITY AND TRIM POOR QUALITY READS ###

# ### change this to make a folder for the chop data

# ### Use porechop to assess read quality and trim the poor quality reads
# ### requires input dir, output file and flag to discard reads with adapters in the middle (this is required for nanopolish)

# for folder in /home3/dd87a/"$project"/demulti_trim/*/ ;

# do

#        folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

#        /software/Porechop-v0.2.4/porechop-runner.py -i "$folder" -o /home3/dd87a/"$project"/demulti_trim/"$folderName"/"$folderName"Chop.fastq --discard_middle

# done


# # ##### ### RUN QC REPORT FOR EACH BARCODE POST PORECHOP ###

# for sample in /home3/dd87a/"$project"/demulti_trim/*/*Chop*.fastq ; 
# #for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

# do

#         folderName=`echo "$sample" | sed -e 's/^.*\(barcode*\/\)\///'| sed -e 's/\.fastq$//'`
        
#         NanoStat --fastq "$sample" --outdir /home3/dd87a/"$project"/statreport -n "$folderName".txt

# done

# ##### POSSIBLY ADD A CHOPPER SECTION IN HERE TO FILTER THE READS TO LONGER READS (>300bp) ???

# ##### ### GENERATE DRAFT GENOME ASSEMBLIES ###

# ### FLYE ASSEMBLY TO USE WITH MEDAKA

# for sample in /home3/dd87a/"$project"/demulti_trim/*/*Chop*.fastq ; 

# do

# 		sample=`echo "$sample" | sed -e 's/^.*\(barcode*\/\)\///' | sed -e 's/\.fastq$//'`

#         sampleID=`echo "$sample" | sed -e 's/^.*trim\/barcode..\///' | sed -e 's/\.fastq$//'`
       
#         /software/flye-v2.9.1/bin/flye --nano-hq "$sample".fastq --out-dir /home3/dd87a/"$project"/draftAssembly/"$sampleID" --genome-size=135000 --iterations 4

# done

# ##### ### ASSESS QUALITY OF THE DE NOVO CONTIGS GENERATED WITH FLYE ###

# ### NOT WORKING QUAST NEEDS PYTHON 2!!

# # for draftAssembly in /home3/dd87a/"$project"/draftAssembly/*/assembly.fasta ; 

# # do

# # 		#sample=`echo "$sample" | sed -e 's/^.*\(barcode*\/\)\///' | sed -e 's/\.fastq$//'`

# # 		assemblyID=`echo "$draftAssembly" | sed -e 's/^.*Assembly\///' | sed -e 's/\/assembly.fastq$//'`
# # 	 # /software/quast-v5.0.2/quast.py <input files> -t 8 -o <outdir> -r <refgenome> <files_with_contigs>

# # 		/software/quast-v5.0.2/quast.py -t 8 -o /home3/dd87a/"$project"/draftAssembly/"$assemblyID"/ -r "$refgenome" 

# # done

# ##### ### FILTER READS TO MAP ###

# ##### reduce the number of really short reads and filter on quality >10

# # for folder in /home3/dd87a/"$project"/demulti_trim/*/;

# # do

# #        folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

# #        Nanofilt 

# # done


# ##### ### MAP READS TO REFERENCE ###

# ##### use minimap2 to map reads to reference genome

# # -ax CIGAR and output to SAM, (x) map-ont - map ONT data, align noisy long reads 10% error rate to ref
# # -bS output to BAM, compatibitity check          -F exclude FLAG'd alignments mapping reads
# # -o output file 

# for folder in /home3/dd87a/"$project"/demulti_trim/*/;

# do

#        folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

#        NanoFilt -l 500 -q 10 /home3/dd87a/"$project"/demulti_trim/"$folderName"/"$folderName"Chop.fastq | minimap2 -ax map-ont "$project"/"$refseq" | samtools view -bS -F 4 | samtools sort -o /home3/dd87a/"$project"/assemblies/"$folderName"draft.bam
#        #minimap2 -ax map-ont "$project"/"$refseq" /home3/dd87a/"$project"/demulti_trim/"$folderName"/"$folderName"Chop.fastq | samtools view -bS -F 4 | samtools sort -o /home3/dd87a/"$project"/assemblies/"$folderName"draft.bam

#        samtools view -h /home3/dd87a/"$project"/assemblies/"$folderName"draft.bam > /home3/dd87a/"$project"/assemblies/"$folderName"draft.sam 

# done



### RUN QC REPORT FOR EACH BARCODE POST ASSEMBLY

# for algn in /home3/dd87a/"$project"/assemblies/*.bam ; 

# do

#         algnName=`echo "$algn" | sed -e 's/\.bam$//'`
        
#         NanoStat --bam "$algn" --outdir /home3/dd87a/"$project"/statreport -n "$algnName".txt

# done


##### ### POLISH WITH MEDAKA AND CREATE CONSENSUS SEQUENCE ###
##### CVR thread limit is 8!

### Porechopped reads : /home3/dd87a/"$project"/demulti_trim/"$folderName"/"$folderName"Chop.fastq
### Reference assembly : /home3/dd87a/"$project"/assemblies/"$folderName"draft.sam 
### Draft assembly contigs : /home3/dd87a/"$project"/draftAssembly/"$sampleID"/ambly.fasta

### Model type : {pore}_{device}_{caller variant}_{caller version}

mkdir /home3/dd87a/"$project"/medaka

for folder in /home3/dd87a/"$project"/demulti_trim/*/;
# /home3/dd87a/"$project"/demulti_trim/*/*Chop*.fastq
do

       folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`
       sampleID=`echo "$folder"/*Chop*.fastq | sed -e 's/^.*trim\/barcode...\///' | sed -e 's/\.fastq$//'`

	   medaka_consensus -i /home3/dd87a/"$project"/demulti_trim/"$folderName"/"$sampleID".fastq -d /home3/dd87a/"$project"/draftAssembly/"$sampleID"/assembly.fasta -o /home3/dd87a/"$project"/medaka/ -t 8 -m r941_min_high_g360

done



} > medakaPipelineLog.txt






## MINIASM PIPELINE WITHOUT FLYE DE NOVO ASSEMBLY ##

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

refLabel=`echo "$refseq" | sed -e 's/\.fa//'`

##### ACTIVATE "ont" CONDA ENV - NEEDED FOR DATA PROCESSING STEPS
source /home3/dd87a/miniconda3/bin/activate ont

mkdir /home3/dd87a/"$project"/basecalled;
mkdir /home3/dd87a/"$project"/demulti_trim;
mkdir /home3/dd87a/"$project"/assemblies;
mkdir /home3/dd87a/"$project"/draftAssembly;
mkdir /home3/dd87a/"$project"/statreport/;
mkdir /home3/dd87a/"$project"/porechop/;
mkdir /home3/dd87a/"$project"/miniasm/;

dateStamp=date +%d%m%y


##### ### UNZIP FAST5 FILES ###

# for arxv in "$project"/allFast5/arxvd/\*.zip;

# 	do
# 		unzip "$arxv" -d $project"/allFast5/

# 	done


##### ### BASECALLING ###

##### BASECALL ON GPU

# /software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/"$project"/allFast5 --save_path /home3/dd87a/"$project"/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96 -x "cuda:0"

# ##### ### DEMULTIPLEXING ###

# ##### DEMULTIPLEX BASECALLED READS

guppy_barcoder -r --input_path /home3/dd87a/"$project"/basecalled/pass --save_path /home3/dd87a/"$project"/demulti_trim --trim_barcodes 

cp -R /home3/dd87a/"$project"/demulti_trim /home3/dd87a/"$project"/demulti_trim_backup

# ##### ### INITIAL QC STEP FOR EACH BARCODE ###

# ##### RUN QC REPORT FOR EACH BARCODE - WORKING 11/12/23

 for folder in /home3/dd87a/"$project"/demulti_trim/*/ ;
### for folder in /home3/dd87a/FHV_iso_26may23/demulti_trim/*/ ;

 do

    folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

# Check this file reference with the asterisk

    NanoStat --fastq "$folder"*.fastq --outdir /home3/dd87a/"$project"/statreport/ -n "$folderName"preChop.txt

 done

##### ### ASSESS READ QUALITY AND TRIM POOR QUALITY READS ###

### Use porechop to assess read quality and trim the poor quality reads
### requires input dir, output file and flag to discard reads with adapters in the middle (this is required for nanopolish)

for folder in /home3/dd87a/"$project"/demulti_trim/*/ ;

do

    folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`

    /software/Porechop-v0.2.4/porechop-runner.py -i "$folder" -o /home3/dd87a/"$project"/porechop/"$folderName"Chop.fastq --discard_middle

done

# ##### ### RUN QC REPORT FOR EACH BARCODE POST PORECHOP ###

for sample in /home3/dd87a/"$project"/porechop/*Chop.fastq ; 

do

    folderName=`echo "$sample" | sed -e 's/^.*\(barcode*\/\)\///'| sed -e 's/\.fastq$//'`
        
    NanoStat --fastq "$sample" --outdir /home3/dd87a/"$project"/statreport/ -n "$folderName"postChop.txt

done

##### ### FILTER READS AND MAP READS TO REFERENCE ->   MAPPED READS TO FOLDER FOR DE NOVO DRAFT ASSEMBLY ( Minimap used as a filter, not to create a mapped assembly ) ###

##### use minimap2 to map reads to reference genome

# -ax CIGAR and output to SAM, (x) map-ont - map ONT data, align noisy long reads 10% error rate to ref
# -bS output to BAM, compatibitity check          -F exclude FLAG'd alignments mapping reads
# -o output file 

for sample in /home3/dd87a/"$project"/porechop/;

do

    folderName=`echo "$folder" | sed -e 's/^.*porechop\///'| sed -e 's/Chop\.fastq$//'`
# Filter reads by length and quality

    NanoFilt -l 300 -q 10 /home3/dd87a/"$project"/porechop/"$folderName"Chop.fastq > /home3/dd87a/"$project"/porechop/"$folderName"Filt.fastq

# Map reads to ref and move mapped reads to file
    
    minimap2 -ax map-ont "$refseq" /home3/dd87a/"$project"/porechop/"$folderName"Filt.fastq | /software/samtools-v1.16.1/bin/samtools view -b -F 4 > /home3/dd87a/"$project"/porechop/"$folderName"_mapped.sam

# Map reads to ref and move unmapped reads to file
    
    minimap2 -ax map-ont "$refseq" /home3/dd87a/"$project"/porechop/"$folderName"Filt.fastq | /software/samtools-v1.16.1/bin/samtools view -b -f 4 > /home3/dd87a/"$project"/porechop/"$folderName"_unmapped.sam

# Convert sam file of mapped reads to fastq to be used as input for de novo assembly
    
    /software/bedtools-v2.29.1/bin/bamToFastq -i /home3/dd87a/"$project"/porechop/"$folderName"_mapped.sam -fq /home3/dd87a/"$project"/porechop/"$folderName"_mapped.fastq 

# Filter mapped reads to minimum reads length - already filtered above using nanofilt

    #/software/bbmap-v38.90/reformat.sh in=/home3/dd87a/"$project"/porechop/"$folderName"_mapped.fastq  out=/home3/dd87a/"$project"/demulti_trim/"$folderName"/"$folderName"_mappedF.fastq  minlength=200

done  

    
for filename in /home3/dd87a/"$project"/porechop/*_mapped.fastq  ; 

    do

        sample= `echo "$filename" | sed -e 's/^.*porechop\///'| sed -e 's/_mapped.fastq$//'
        sampleShort=`echo "$sample" | sed -e 's/barcode/bc/'

        mkdir /home3/dd87a/"$project"/miniasm/"$sample"

        minimap2 -x ava-ont "$filename" "$filename" | gzip -1 > /home3/dd87a/"$project"/miniasm/"$sample"/minimap.paf.gz

        # build "unitigs" - (high confidence contigs)

        miniasm -f "$filename" /home3/dd87a/"$project"/miniasm/"$sample"/minimap.paf.gz > /home3/dd87a/"$project"/miniasm/"$sample"/miniasm.gfa

        # convert to fasta

        awk '/^S/{print ">"$2"\n"$3}' /home3/dd87a/"$project"/miniasm/"$sample"/miniasm.gfa > /home3/dd87a/"$project"/miniasm/"$sample"/miniasm.fasta

        # generate polished consensus

        mkdir /home3/dd87a/"$project"/medaka/"$sample"/miniasm_consensus_"$sampleShort"

        medaka_consensus -i "$filename" -d /home3/dd87a/"$project"/miniasm/"$sample"/miniasm.fasta -o /home3/dd87a/"$project"/medaka/miniasm_consensus_"$sampleShort" -t 8 -m r941_min_high_g360

        # there are two unitigs after this step!

        #/software/flye-v2.9.1/bin/flye --nano-hq /home3/dd87a/fhv_02aug23/medaka/miniasm_consensus_bc01/miniasm_consensus_bc01.fasta --out-dir /home3/dd87a/fhv_02aug23/draftAssembly/barcode01_miniCons --genome-size=135000 --iterations 4

        # full genome assembly asm5/10/20 long read assembly to ref presets - with divergence accounted for 5/10/20% 

        minimap2 -ax asm5 "$refseq" /home3/dd87a/"$project"/medaka/miniasm_consensus_"$sampleShort"/miniasm_consensus_"$sampleShort".fasta | /software/samtools-v1.16.1/bin/samtools view -bS -F 4 | /software/samtools-v1.16.1/bin/samtools sort -o /home3/dd87a/"$project"/assemblies/miniasmwg"$sampleShort".bam

        # index bam file

        /software/samtools-v1.16.1/bin/samtools index /home3/dd87a/"$project"/assemblies/miniasmwg"$sampleShort".bam

        # Need to generate a consensus from the bam files 

        # alfred consensus -t ill -f bam -p chr4:500500 input.bam
        # bcftools consensus ?

    done








} > miniasmPipelineLog_"$dateStamp".txt
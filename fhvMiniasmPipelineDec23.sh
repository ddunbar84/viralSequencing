
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

##### ### FILTER READS AND MAP READS TO REFERENCE ->   MAPPED READS TO FOLDER FOR DE NOVO DRAFT ASSEMBLY ###

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

## Filter mapped reads to minimum reads length - already filtered above using nanofilt

    #/software/bbmap-v38.90/reformat.sh in=/home3/dd87a/"$project"/porechop/"$folderName"_mapped.fastq  out=/home3/dd87a/"$project"/demulti_trim/"$folderName"/"$folderName"_mappedF.fastq  minlength=200

done    

##### ### GENERATE DRAFT DE NOVO GENOME ASSEMBLIES AS INPUT FOR MEDAKA ###

### FLYE ASSEMBLY TO USE WITH MEDAKA

#for sample in /home3/dd87a/"$project"/demulti_trim/*/*_mappedF.fastq ; 
for sample in /home3/dd87a/"$project"/porechop/*.fastq ; 

do

	sample=`echo "$sample" | sed -e 's/^.*porechop\///' | sed -e 's/\.fastq$//'`

    sampleID=`echo "$sample" | sed -e 's/^.*\/barcode..\///' | sed -e 's/_mapped\.fastq$//'`
       
    /software/flye-v2.9.1/bin/flye --nano-hq "$sample"_mapped.fastq --out-dir /home3/dd87a/"$project"/draftAssembly/"$sampleID" --genome-size=135000 --iterations 4

done

# Rename files and move folder, so they are all in one folder location for next steps

for assembly in /home3/dd87a/"$project"/draftAssembly/*/assembly.fasta;

do 

	sampleID=`echo "$assembly" | sed -e 's/^.*Assembly\///' | sed -e 's/Chop.*.fasta$//'`

	cp "$assembly" /home3/dd87a/"$project"/draftAssembly/"$sampleID"Assembly.fastq

done


##### ### ASSESS QUALITY OF THE DE NOVO CONTIGS GENERATED WITH FLYE ###

### NOT WORKING QUAST NEEDS PYTHON 2!!

# for draftAssembly in /home3/dd87a/"$project"/draftAssembly/*/assembly.fasta ; 

# do

# 		#sample=`echo "$sample" | sed -e 's/^.*\(barcode*\/\)\///' | sed -e 's/\.fastq$//'`

# 		assemblyID=`echo "$draftAssembly" | sed -e 's/^.*Assembly\///' | sed -e 's/\/assembly.fastq$//'`
# 	 # /software/quast-v5.0.2/quast.py <input files> -t 8 -o <outdir> -r <refgenome> <files_with_contigs>

# 		/software/quast-v5.0.2/quast.py -t 8 -o /home3/dd87a/"$project"/draftAssembly/"$assemblyID"/ -r "$refgenome" 

# done


##### ### ALIGN DE NOVO ASSEMBLIES TO REFERENCE TO CHECK QUALITY OF THE DE NOVO ###assembly

##### Individually

# for assembly in /home3/dd87a/"$project"/draftAssembly/*/assembly.fasta;

# do 

# 	sampleID=`echo "$assembly" | sed -e 's/^.*Assembly\///' | sed -e 's/Chop.*.fasta$//'`
# 	cat "$assembly" /home3/dd87a/NC_013590.fa > /home3/dd87a/"$project"/draftAssembly/"$sampleID"Merged.fasta
# 	mafft /home3/dd87a/"$project"/draftAssembly/"$sampleID"Merged.fasta > /home3/dd87a/"$project"/draftAssembly/"$sampleID"Aligned.fasta

# done

##### Align all

for assembly in /home3/dd87a/"$project"/draftAssembly/*Assembly.fastq;

# Need to label contigs in the assemblies otherwise Medaka throws an error

do 
	sampleID=`echo "$assembly" | sed -e 's/^.*Assembly\///' | sed -e 's/Assembly.*.fastq$//'`
	sed "s/>contig/>${sampleID}contig/" "$assembly" > /home3/dd87a/"$project"/draftAssembly/"$sampleID"LabAsmbl.fastq 

done

cat /home3/dd87a/"$project"/draftAssembly/*LabAsmbl.fastq /home3/dd87a/NC_013590.fa > /home3/dd87a/"$project"/draftAssembly/allsamplesMerged.fasta
mafft /home3/dd87a/"$project"/draftAssembly/allsamplesMerged.fasta > /home3/dd87a/"$project"/draftAssembly/allsamplesAligned.fasta


### RUN QC REPORT FOR EACH BARCODE POST ASSEMBLY

# for algn in /home3/dd87a/"$project"/assemblies/*.bam ; 

# do

#         algnName=`echo "$algn" | sed -e 's/\.bam$//'`
#         NanoStat --bam "$algn" --outdir /home3/dd87a/"$project"/statreport -n "$algnName".txt        


# done


##### POLISH WITH MEDAKA AND CREATE CONSENSUS SEQUENCE #####

mkdir /home3/dd87a/"$project"/medaka

for folder in /home3/dd87a/"$project"/porechop/*/;
# /home3/dd87a/"$project"/porechop/*/*Chop*.fastq
do

    folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`
    sampleID=`echo "$folder"/*Chop*.fastq | sed -e 's/^.*trim\/barcode...\///' | sed -e 's/\.fastq$//'`

    mkdir /home3/dd87a/"$project"/medaka/"$sampleID";

	medaka_consensus -i /home3/dd87a/"$project"/porechop/"$folderName"/"$sampleID".fastq -d /home3/dd87a/"$project"/draftAssembly/"$sampleID"/assembly.fasta -o /home3/dd87a/"$project"/medaka/"$sampleID" -t 8 -m r941_min_high_g360


done

for consensus in /home3/dd87a/"$project"/medaka/consensus/*consensus.fasta;

do 
        sampleID=`echo "$consensus" | sed -e 's/^.*consensus\///' | sed -e 's/consensus.fasta$//'`
        echo "$sampleID"
# change "contig" label to "sampleID_contig"    
        sed "s/>contig/>${sampleID}contig/" "$consensus" > /home3/dd87a/"$project"/medaka/"$sampleID"LabCons.fasta 

done


# All merged and labelled, then rename the refseq without punctuation
cat /home3/dd87a/"$project"/medaka/*.fasta "$refseq" > /home3/dd87a/"$project"/medaka/allMergedLab.fasta
sed 's/>ref.*/>NC013590/1' /home3/dd87a/"$project"/medaka/allMergedLab.fasta > /home3/dd87a/"$project"/medaka/allMergedLab_edit.fasta 



### Rename concatenated consensus
# sed 's/>ref.*/>NC013590/1' /home3/dd87a/fhv_02aug23/assemblies/allMerged.fasta > /home3/dd87a/fhv_02aug23/assemblies/allMerged_edit.fasta 
# mafft /home3/dd87a/fhv_02aug23/medaka/consensus/ > /home3/dd87a/fhv_02aug23/medaka/allAligned.fasta

##### QUALITY ASSESS WITH MEDAKA ???? ######




##### VARIANT ANALYSIS USING MEDAKA VARIANT CALLER  ######

for folder in /home3/dd87a/"$project"/porechop/*/;
# /home3/dd87a/"$project"/porechop/*/*Chop*.fastq
do

    folderName=`echo "$folder" | sed -e 's/^.*\(trim\)\///'| sed -e 's/\/$//'`
    sampleID=`echo "$folder"/*Chop*.fastq | sed -e 's/^.*trim\/barcode...\///' | sed -e 's/\.fastq$//'`

    medaka_haploid_variant -i /home3/dd87a/"$project"/porechop/"$folderName"/"$sampleID".fastq -r "$project"/"$refseq" -o /home3/dd87a/"$project"/medaka/"$sampleID"VC

done

##### RENAME ALL MEDAKA OUTPUT FILES CONSENSUS AND VC  ######

##### CONSENSUS  ######

for filename in /home3/dd87a/"$project"/medaka/*Chop/* ; 

	do

		sampleID=`echo "$filename" | sed -e 's/^.*\(barcode\)/bc/'| sed -e 's/Chop.*\/.*$//'`
		filepath=`echo "$filename" | sed -e "s/Chop\/.*$/Chop\//"`
		filename2=`echo "$filename" | sed -e "s/^.*Chop\///"`

		mv "$filename" "$filepath"/"$sampleID"_"$filename2";

	done

##### VC  ######

for filename in /home3/dd87a/"$project"/medaka/*ChopVC/* ; 

	do

		sampleID=`echo "$filename" | sed -e 's/^.*\(barcode\)/bc/'| sed -e 's/Chop.*\/.*$//'`
		filepath=`echo "$filename" | sed -e "s/ChopVC\/.*$/ChopVC\//"`
		filename2=`echo "$filename" | sed -e "s/^.*ChopVC\///"`

		mv "$filename" "$filepath"/"$sampleID"_"$filename2";

	done		

##### ### MINIASM PIPELINE FOR ASSEMBLY ###

# map reads to themselves

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



##### MAXIMUM LIKELIHOOD TREE OF CONSENSUS GENOMES  ######

### needs to have all consensus seqs and ref genome concatenated into a single FASTA file.


#/software/raxml-ng-1.0.2/bin/raxml-ng-mpi --search1 --msa fhvGenbank_alignStripped --model GTR+G

} > medakaPipelineLog_"$dateStamp".txt





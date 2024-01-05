
### !!!!! need to add a Nanofilt step - but this has already been done with the files I'm working with! 

# Map reads to ref and move mapped reads to file

minimap2 -ax map-ont fhv_02aug23/NC_013590.fa /home3/dd87a/fhv_02aug23/demulti_trim/barcode01/barcode01Filt.fastq | /software/samtools-v1.16.1/bin/samtools view -b -F 4 > /home3/dd87a/fhv_02aug23/demulti_trim/barcode01/barcode01_mapped.sam

# Map reads to ref and move unmapped reads to file

minimap2 -ax map-ont fhv_02aug23/NC_013590.fa /home3/dd87a/fhv_02aug23/demulti_trim/barcode01/barcode01Filt.fastq | /software/samtools-v1.16.1/bin/samtools view -b -f 4 > /home3/dd87a/fhv_02aug23/demulti_trim/barcode01/barcode01_unmapped.sam

# Convert sam file of mapped reads to fastq to be used as input for de novo assembly

/software/bedtools-v2.29.1/bin/bamToFastq -i fhv_02aug23/demulti_trim/barcode01/barcode01_mapped.sam -fq fhv_02aug23/demulti_trim/barcode01/barcode01_mapped.fastq 

# Filter minimum read length

/software/bbmap-v38.90/reformat.sh in=fhv_02aug23/demulti_trim/barcode01/barcode01_mapped.fastq  out=fhv_02aug23/demulti_trim/barcode01/barcode01_mappedF.fastq  minlength=200

### For script - need to concat files and reference and then align - if looking at the ref mapped assembly AND have run a samtools consensus!!!

/software/flye-v2.9.1/bin/flye --nano-hq fhv_02aug23/demulti_trim/barcode01/barcode01_mapped.fastq --out-dir /home3/dd87a/fhv_02aug23/draftAssembly/barcode01_mapped --genome-size=135000 --iterations 4

# /software/flye-v2.9.1/bin/flye --nano-hq fhv_02aug23/demulti_trim/barcode01/barcode01_mapped_Copy.fastq --out-dir /home3/dd87a/fhv_02aug23/draftAssembly/barcode01_mapped --genome-size=135000 --iterations 2 --scaffold   --threads 8


# map reads to themselves

minimap2 -x ava-ont fhv_02aug23/demulti_trim/barcode01/barcode01_mapped_Copy.fastq fhv_02aug23/demulti_trim/barcode01/barcode01_mapped_Copy.fastq | gzip -1 > fhv_02aug23/demulti_trim/barcode01/minimap.paf.gz

# build "unitigs" - (high confidence contigs)

miniasm -f fhv_02aug23/demulti_trim/barcode01/barcode01_mapped_Copy.fastq fhv_02aug23/demulti_trim/barcode01/minimap.paf.gz > fhv_02aug23/demulti_trim/barcode01/miniasm.gfa

# convert to fasta

awk '/^S/{print ">"$2"\n"$3}' fhv_02aug23/demulti_trim/barcode01/miniasm.gfa > fhv_02aug23/demulti_trim/barcode01/miniasm.fasta

# generate polished consensus

medaka_consensus -i fhv_02aug23/demulti_trim/barcode01/barcode01_mapped_Copy.fastq -d fhv_02aug23/demulti_trim/barcode01/miniasm.fasta -o /home3/dd87a/fhv_02aug23/medaka/miniasm_consensus_bc01 -t 8 -m r941_min_high_g360

# there are two unitigs after this step!

#/software/flye-v2.9.1/bin/flye --nano-hq /home3/dd87a/fhv_02aug23/medaka/miniasm_consensus_bc01/miniasm_consensus_bc01.fasta --out-dir /home3/dd87a/fhv_02aug23/draftAssembly/barcode01_miniCons --genome-size=135000 --iterations 4

# full genome assembly asm5/10/20 long read assembly to ref presets - with divergence accounted for 5/10/20% 

minimap2 -ax asm5 fhv_02aug23/NC_013590.fa /home3/dd87a/fhv_02aug23/medaka/miniasm_consensus_bc01/miniasm_consensus_bc01.fasta | /software/samtools-v1.16.1/bin/samtools view -bS -F 4 | /software/samtools-v1.16.1/bin/samtools sort -o fhv_02aug23/assemblies/miniasmwgbc01.bam

# index bam file

/software/samtools-v1.16.1/bin/samtools index fhv_02aug23/assemblies/miniasmwgbc01asm10.bam

# call consensus

/software/samtools-v1.16.1/bin/samtools consensus -f fastq fhv_02aug23/assemblies/miniasmwgbc01asm10.bam -o fhv_02aug23/assemblies/miniasmCons_bc01asm10.fastq

# /software/samtools-v1.16.1/bin/samtools consensus -f fastq /home3/dd87a/fhv_02aug23/assemblies/"$folderName"RA.bam -o /home3/dd87a/fhv_02aug23/assemblies/miniasm_fullCons.fastq




minipolish -t 2 nanofilt_result.fastq miniasm.gfa > minipolished_assembly.gfa
dnadiff -p dnadiff ~/course_data/precompiled/chr17.fasta miniasm.fasta


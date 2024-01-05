#!bin/bash 

#tmux
#run process
#ctrl +b, d (detach from tmux)

#in shell tmux attach to get back into tmux screen

#Log in to Alpha2
#screen -S GPU3
#ssh gpu3 (prompts for host then password)
#nvidia-smi

# use for RNA PCR barcodes
guppy_basecaller --input_path /home3/dd87a/fast5folder/ --save_path /home3/dd87a/fast5folder/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109
# basecall on GPU
/software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/fast5folder/ --save_path /home3/dd87a/fast5folder/basecalled --flowcell FLO-MIN106 --kit SQK-PCB109 -x "cuda:0"

# use for DNA barcodes
guppy_basecaller --input_path /home3/dd87a/fast5folder/ --save_path /home3/dd87a/fast5folder/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96
# basecall on GPU
/software/ont-guppy-GPU-v4.5.3/bin/guppy_basecaller --input_path /home3/dd87a/fast5folder/ --save_path /home3/dd87a/fast5folder/basecalled --flowcell FLO-MIN106 --kit SQK-RBK110-96 -x "cuda:0"

# demultiplex reads
guppy_barcoder -r --input_path /home3/dd87a/fast5folder/basecalled_pass --save_path /home3/dd87a/fast5folder/basecalled/demultiplexed_trimmed --trim_barcodes 
#guppy_barcoder -r --input_path /home3/dd87a/FHV_test_28feb23/basecalled/pass --save_path /home3/dd87a/FHV_test_28feb23/basecalled/pass/demultiplexed_trimmed --trim_barcodes 

# Add in a QC check step?

# Use porechop to assess read quality and trim the poor quality reads
# input, output and discard reads with adapters in the middle (required for nanopolish)
/software/Porechop-v0.2.4/porechop-runner.py -i /home3/dd87a/fast5folder/basecalled/demultiplexed_trimmed -o qtrimmedReads.fastq --discard_middle

# Add in a QC check step?

# use minimap2 (or flye) to map reads to reference genome
# -ax CIGAR and output to SAM, (x) map-ont - map ONT data, align noisy long reads 10% error rate to ref
# -bS output to BAM, compatibitity check          -F exclude FLAG'd alignments mapping reads
# -o output file 
minimap2 -ax map-ont reference.fa runA.fastq | samtools view -bS -F 4 | samtools sort -o runAtoRef.bam

# Add a WIMP mapping protocol?

# Add in a QC check step?

nanopolish-v0.14.0  

/software/nanopolish-v0.14.0/nanopolish

nanopolish vcf2fasta -g draft.fa polished.*.vcf > polished_genome.fa

# ivar call consensus / variant calling?



#!bin/bash

##############
##Dmitry Sutormin, 2024##
##General read mapping and coverage depth preparation##

#Shell script that makes QC of the reads before and after the trimming procedure. 
#Than script maps trimmed only paired reads to the reference genome, prepares
#sorted and indexed BAM-files suitable for visualization with IGV

#Requirements: factqc, bwa mem, samtools 
##############


#######
#Variables to be defined.
#######

#Path to the working directory, contains /Raw_data folder with raw reads files in fq.gz format.
PWD='/local/DSutormin/rRNA_depletion_project/PAO1'
echo $PWD
cd $PWD

#Path to the reference genome in fasta/fna format.
Ref_genome='/local/DSutormin/rRNA_depletion_project/PAO1/Reference_genome/Pseudomonas_aeruginosa_pao1_GCA_000006765.1_ASM676v1_genomic.fna'
#Path to fastqc executable.
fastqc='/local/DSutormin/Programms/FastQC/fastqc'
#Path to bwa executable.
bwa='/local/DSutormin/Programms/bwa'
#Path to samtools executable.
samtools='/local/DSutormin/Programms/samtools'


#Initial quality control
echo '
#######################
Initial quality control is in progress...
#######################
'
mkdir $PWD/Fastqc_analysis/
$fastqc -t 20 -o $PWD/Fastqc_analysis/ $PWD/Raw_data/*


#Reads mapping to the reference genome: make SAM-files
echo '
#######################
Reads mapping, SAM files generation...
#######################
'
mkdir $PWD/SAM/
$bwa index $Ref_genome
for i in `ls -a $PWD/Raw_data/ | grep 'fq.gz' | sed -r "s/(.+)_[1,2]\.fq\.gz/\1/g" | uniq | sort -d`; do 
echo ${i}
$bwa mem -t 20 $Ref_genome $PWD/Raw_data/${i}_1.fq.gz > $PWD/SAM/$i.sam; done


#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
echo '
#######################
BAM files preparation...
#######################
'
mkdir $PWD/BAM/
for i in `ls -a $PWD/SAM/ | grep '.sam' | sed -r "s/(.+).sam/\1/g"`; do 
$samtools view -@ 10 -S -b $PWD/SAM/${i}.sam > $PWD/BAM/${i}.bam ; done


#Sorts BAM-files
echo '
#######################
BAM files sorting...
#######################
'
mkdir $PWD/BAM_sorted/
for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
$samtools sort -@ 10 $PWD/BAM/${i}.bam -o $PWD/BAM_sorted/${i}_sorted.bam ; done


#Makes indexes for bam files.
echo '
#######################
BAM files indexing...
#######################
'
for i in `ls -a $PWD/BAM_sorted/`; do 
$samtools index -@ 10 $PWD/BAM_sorted/${i} ; done


#Creates coverage depth file (bed-like format).
echo '
#######################
Coverage depth calculation...
#######################
'
mkdir $PWD/Genome_coverage_depth/
for i in `ls -a $PWD/BAM_sorted/ | grep '.bam$' | sed -r "s/(.+).bam/\1/g"`; do 
$samtools depth -@ 10 -a -o $PWD/Genome_coverage_depth/${i}.bed $PWD/BAM_sorted/${i}.bam ; done

echo '
Script finished succesfully!'

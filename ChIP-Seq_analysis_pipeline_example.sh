#!bin/bash

##############
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Shell script that downloads data from ncbi, performs QC of the reads before and after the trimming procedure. 
#Than script maps trimmed only paired reads to the reference genome, prepares
#sorted and indexed BAM-files suitable for visualization with IGV

#Requirements: sra toolkit, factqc, trimmomatic, bwa mem, samtools 
#This variables should be in the path (or replace them with the path to the particular program)
##############



#######
#Variables to be defined.
#######

#Path to the working directory, contains /Raw_data folder with raw reads files.
PWD='/home/cls01/Data_Hi-C/Topo-Seq_data/TopoI_ChIP-Seq/Ec_TopoI'
echo $PWD
cd $PWD

#Path to SRA folder.
#SRA_path='/home/cls01/ncbi/public/sra/'

#Path to the file containing sequencing adapters sequences for trimmomatic uses. Typically in the Trimmomatic-0.36/adapters/XXX.fa
Adapters='/home/cls01/Prog/Trimmomatic-0.38/adapters/All_TruSeq.fa'
trimmomatic='/home/cls01/Prog/Trimmomatic-0.38/trimmomatic-0.38.jar'
#Path to the reference genome.
#Ref_genome_link='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/767/705/GCA_000767705.1_ASM76770v1/GCA_000767705.1_ASM76770v1_genomic.fna.gz'
Ref_genome=$PWD/Genome/E_coli_w3110_G_Mu.fasta



#######
#Download data from ncbi and prepare fastq files.
#######
#echo '
########################
#NGS data fetching from NCBI, fastq preparation...
########################
#'

#Clean the SRA directory.
#rm $SRA_path/SRR*

#Download data from NCBI
#prefetch SRR8149481 SRR8149480 SRR8149482 SRR8149479

#Prepare fastq files.
#mkdir $PWD/Raw_data
#fastq-dump --split-files --outdir $PWD/Raw_data $SRA_path/SRR*


#######
#Quality control and sequencing data preparation.
#######
echo '
#######################
Initial quality control is in progress...
#######################
'

#Initial quality control
mkdir $PWD/Fastqc_analysis/
mkdir $PWD/Fastqc_analysis/Initial
fastqc -t 40 -o $PWD/Fastqc_analysis/Initial $PWD/Raw_data/*



#######
#Reads trimming
#######
echo '
#######################
Reads trimming...
#######################
'

mkdir $PWD/Trimmed/
for i in `ls -a $PWD/Raw_data/ | grep 'fastq' | sed -r "s/(.+)_R[1,2]_001\.fastq\.gz/\1/g" | uniq | sort -d`; do
echo $i
java -jar $trimmomatic PE -threads 40 -phred33 $PWD/Raw_data/${i}_R1_001.fastq.gz $PWD/Raw_data/${i}_R2_001.fastq.gz $PWD/Trimmed/${i}_paired_R1.fastq.gz $PWD/Trimmed/${i}_unpaired_R1.fastq.gz $PWD/Trimmed/${i}_paired_R2.fastq.gz $PWD/Trimmed/${i}_unpaired_R2.fastq.gz ILLUMINACLIP:$Adapters:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:30 ; done



#######
#Quality control after the trimming procedure
#######
echo '
#######################
Quality control after trimming...
#######################
'

mkdir $PWD/Fastqc_analysis/Trimmed/
fastqc -t 40 -o $PWD/Fastqc_analysis/Trimmed/ $PWD/Trimmed/*



#######
#Download and prepare index for reference genome.
#######
#echo '
########################
#Reference genome fetching, unzipping, indexing...
########################
#'

#mkdir $PWD/Ref_genome
#wget -O $PWD/Ref_genome/$Ref_genome.gz $Ref_genome_link
#gzip -d $PWD/Ref_genome/$Ref_genome.gz
#bwa index $PWD/Ref_genome/$Ref_genome


#######
#Reads mapping, alignment conversion to IGV-compatible format (sorted indexed BAM).
#######
echo '
#######################
Reads mapping, SAM files generation...
#######################
'

#Reads mapping to the reference genome: make SAM-files
mkdir $PWD/SAM/
for i in `ls -a $PWD/Trimmed/ | grep '_paired_' | sed -r "s/(.+)_paired_R[1,2]\.fastq\.gz/\1/g" | uniq | sort -d`; do 
bwa mem -t 40 $Ref_genome $PWD/Trimmed/${i}_paired_R1.fastq.gz $PWD/Trimmed/${i}_paired_R2.fastq.gz > $PWD/SAM/$i.sam; done



#######
#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
#######
echo '
#######################
BAM files preparation...
#######################
'
mkdir $PWD/BAM/
#Makes BAM-files
for i in `ls -a $PWD/SAM/ | grep '.sam' | sed -r "s/(.+).sam/\1/g"`; do 
samtools view -S -b $PWD/SAM/${i}.sam > $PWD/BAM/${i}.bam ; done

#Sorts BAM-files
echo '
#######################
BAM files sorting...
#######################
'
mkdir $PWD/BAM_sorted/
for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools sort $PWD/BAM/${i}.bam -o $PWD/BAM_sorted/${i}_sorted.bam ; done


#Converts bam to bed.
echo '
#######################
BAM to bed conversion...
#######################
'

mkdir $PWD/Cov_depth
for i in `ls -a $PWD/BAM_sorted/ | grep '.bam' | sed -r "s/(.+)_sorted.bam/\1/g"`; do
samtools depth -a $PWD/BAM_sorted/${i}_sorted.bam -o $PWD/Cov_depth/${i}.bed; done



#Removes PCR-duplicates.
echo '


#######################
PCR-duplicates removing...
#######################


'


#Sorts BAM-files by name
echo '
#######################
BAM files sorting by name...
#######################
'
mkdir $PWD/BAM_name_sorted
for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools sort -n $PWD/BAM/${i}.bam -o $PWD/BAM_name_sorted/${i}_ns.bam ; done

#Fixmates.
echo '
#######################
Fixmating...
#######################
'
mkdir $PWD/BAM_fixmated
for i in `ls -a $PWD/BAM_name_sorted/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools fixmate -m $PWD/BAM_name_sorted/${i}.bam $PWD/BAM_fixmated/${i}_fm.bam; done

#Sorting by position.
echo '
#######################
BAM sorting by position...
#######################
'
mkdir $PWD/BAM_ns_fm_ps
for i in `ls -a $PWD/BAM_fixmated/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools sort $PWD/BAM_fixmated/${i}.bam -o $PWD/BAM_ns_fm_ps/${i}_ps.bam; done

#Removes duplicates.
echo '
#######################
PCR duplicates removing...
#######################
'
mkdir $PWD/BAM_nodup
for i in `ls -a $PWD/BAM_ns_fm_ps/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do 
samtools markdup $PWD/BAM_ns_fm_ps/${i}.bam $PWD/BAM_nodup/${i}_nodup.bam; done

#Converts bam to bed.
echo '
#######################
BAM to bed conversion...
#######################
'

mkdir $PWD/Cov_depth_nodup
for i in `ls -a $PWD/BAM_nodup/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
samtools depth -a $PWD/BAM_nodup/${i}.bam -o $PWD/Cov_depth_nodup/${i}.bed; done


#Prepare tar.gz archives.
echo '
#######################
tar.gz archives preparation...
#######################
'

tar -czvf Fastqc_analysis_Ec_TopoI.tar.gz $PWD/Fastqc_analysis/
tar -czvf Cov_depth_Ec_TopoI.tar.gz $PWD/Cov_depth/
tar -czvf Cov_depth_nodup_Ec_TopoI.tar.gz $PWD/Cov_depth_nodup/




#Makes index files
echo '
#######################
BAM files indexing...
#######################
'
for i in `ls -a $PWD/BAM_sorted/`; do 
samtools index $PWD/BAM_sorted/${i} ; done

for i in `ls -a $PWD/BAM_nodup/`; do 
samtools index $PWD/BAM_nodup/${i} ; done


#######
#Peak-calling with MACS2
#######
#echo '
########################
#Peak-calling with MACS2...
########################
#'
#Activate virtual environment with python 2.7.9 and MACS2 installed.
#conda activate MACS2_python_2.7.9
#Peak-calling
#mkdir $PWD/Peak_calling
#mkdir $PWD/Peak_calling/FC
#macs2 callpeak -t $PWD/BAM/SRR8149480.bam -c $PWD/BAM/SRR8149481.bam -n TopoI_FC --outdir $PWD/Peak_calling/FC -f BAM -g 4.5e6 --call-summits -B -q 0.01
#mkdir $PWD/Peak_calling/IT
#macs2 callpeak -t $PWD/BAM/SRR8149482.bam -c $PWD/BAM/SRR8149481.bam -n TopoI_IT --outdir $PWD/Peak_calling/FC -f BAM -g 4.5e6 --call-summits -B -q 0.01
#mkdir $PWD/Peak_calling/D108A
#macs2 callpeak -t $PWD/BAM/SRR8149479.bam -c $PWD/BAM/SRR8149481.bam -n TopoI_D108A --outdir $PWD/Peak_calling/FC -f BAM -g 4.5e6 --call-summits -B -q 0.01
#Deactivate virtual environment.
#conda deactivate

echo '
Script ended its work succesfully!
'

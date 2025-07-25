# This document contains all the code used to analyze 2bRAD sequencing data for Jaskiel et al., 2025
# Author: Jacob Jaskiel

# I have adapted much of what is here from Misha Matz's readme file, found here: https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.sh

#================================== Getting Set Up ==========================================
# (Download scripts from Misha Matz's github: https://github.com/z0on/2bRAD_denovo)

# Clone github into working directory:
git clone https://github.com/z0on/2bRAD_denovo.git

# Make everything executable (chmod) -> This way we can access everything. Just copy and paste these commands into the terminal, they will automatically run
chmod +x *.pl
chmod +x *.py
chmod +x *.R
chmod +x *.txt

# Install modules - we will need perl, bowtie2, samtools, and picard. They are pre-installed as modules on TACC and they are already installed on the SCC
# Copy and paste into terminal
module load perl
module load bowtie2
module load samtools
module load picard

# Reference Genome found here: https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_914725855.1/
# SSID fasta file
GCF_914725855.1_fThuAlb1.1_genomic.fasta

# Index transcriptome for bowtie2 mapper (submit as a job as this can take a while)
>indexing
nano indexing 

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea

module load bowtie2
export GENOME_FASTA=GCF_914725855.1_fThuAlb1.1_genomic.fasta
bowtie2-build $GENOME_FASTA $GENOME_FASTA 

qsub indexing


>indexing2
nano indexing2

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing2 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea

module load samtools
samtools faidx $GENOME_FASTA

qsub indexing2


#================================== Demultiplexing and Trimming Barcodes/Poor Quality End Reads ==========================================
# Step 1: Splitting by in-read barcode, deduplicating and quality-filtering the reads

pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/2bRAD_2024

gunzip *.gz

# Creating a file of commands to run (assuming reads are in fastq files, one file per sample.)
/projectnb/mullenl/jaskiel/2bRAD_fastq/2bRAD_2024/2bRAD_trim_launch_dedup_old.pl fastq > trims   #Use the _old script if using oligos purchased before March 2022 like me

# Modify this trims file to be able to submit it as a job on the SCC, designate where the perl script is, and add '&'s to the end of the line
# It looks like this:

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trims # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m be

module load perl

../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_15_S15_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_7_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_14_S14_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_23_S7_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_5_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_16_S16_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_18_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_13_S13_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_8_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_4_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_21_S5_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_11_S11_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_6_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_10_S10_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_20_S4_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_24_S8_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_9_S9_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_1_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_17_S1_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_3_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_12_S12_L002_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_22_S6_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_19_S3_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
../2bRAD_2024/trim2bRAD_2barcodes_dedup_old.pl input=Pool_2_S2_L001_R1_001.fastq site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}" adaptor="AGATC?" sampleID=100 deduplicate=1 bc=[ATGC]{4} &
wait

qsub trims


# Check to see if you have expected number of *.tr0 files
ls -l *.tr0 | wc -l


# Cutadapt info:
# -q is a quality filter, used to trim low-quality ends from reads.
# The comma separated argument here trims the 5' end with a cutoff of 15 and the 3' end with a cutoff of 15.
# -m is a minimum length filter, it discards processed reads that are shorter than LENGTH provided (here, 25).

# For reference-based analysis: trimming poor quality bases off ends:
>trimse

for file in *.tr0; do
echo "cutadapt -q 15,15 -m 25 -o ${file/.tr0/}.trim $file > ${file}_trimlog.txt" >> trimse;  #This is run outside of the trimse job, and fills in trimse with the proper commands but you still need to add the header and load the relevant modules
done

nano trimse

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N trimse # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -j y # Join standard output and error to a single file
#$ -o trimse.qlog # Name the file where to redirect standard output and error

module load cutadapt 
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_ACCA.trim Pool_10_S10_L002_R1_001_ACCA.tr0 > Pool_10_S10_L002_R1_001_ACCA.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_AGAC.trim Pool_10_S10_L002_R1_001_AGAC.tr0 > Pool_10_S10_L002_R1_001_AGAC.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_AGTG.trim Pool_10_S10_L002_R1_001_AGTG.tr0 > Pool_10_S10_L002_R1_001_AGTG.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_CATC.trim Pool_10_S10_L002_R1_001_CATC.tr0 > Pool_10_S10_L002_R1_001_CATC.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_CTAC.trim Pool_10_S10_L002_R1_001_CTAC.tr0 > Pool_10_S10_L002_R1_001_CTAC.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_GACT.trim Pool_10_S10_L002_R1_001_GACT.tr0 > Pool_10_S10_L002_R1_001_GACT.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_GCTT.trim Pool_10_S10_L002_R1_001_GCTT.tr0 > Pool_10_S10_L002_R1_001_GCTT.tr0_trimlog.txt
cutadapt -q 15,15 -m 25 -o Pool_10_S10_L002_R1_001_GTGA.trim Pool_10_S10_L002_R1_001_GTGA.tr0 > Pool_10_S10_L002_R1_001_GTGA.tr0_trimlog.txt
#only first 6 rows shown

# Submit the trimse job:
qsub trimse

# Do we have expected number of *.trim files created?
ls -l *.trim | wc -l #yes, proceed!


#================================== Mapping to Reference Genome ==========================================

export GENOME_FASTA=GCF_914725855.1_fThuAlb1.1_genomic.fasta

# Create a file called 'maps'
>maps
../2bRAD_2024/2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > maps

# Same as you did above for trimse, nano into maps once you created it and add your header and your 'module load line'.

nano maps

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N maps # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o maps.qlog # Name the file where to redirect standard output and error

module load bowtie2
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_17_S1_L001_R1_001_AGAC.trim -S Pool_17_S1_L001_R1_001_AGAC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_11_S11_L002_R1_001_TGTC.trim -S Pool_11_S11_L002_R1_001_TGTC.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_17_S1_L001_R1_001_GCTT.trim -S Pool_17_S1_L001_R1_001_GCTT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_15_S15_L002_R1_001_TGGT.trim -S Pool_15_S15_L002_R1_001_TGGT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_5_S5_L001_R1_001_GACT.trim -S Pool_5_S5_L001_R1_001_GACT.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_21_S5_L001_R1_001_ACCA.trim -S Pool_21_S5_L001_R1_001_ACCA.trim.bt2.sam
bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x GCF_914725855.1_fThuAlb1.1_genomic.fasta -U Pool_14_S14_L002_R1_001_TCAC.trim -S Pool_14_S14_L002_R1_001_TCAC.trim.bt2.sam


# Setting minimum score that a read has to achieve if the read is 33 base pairs in length (check this) so minimum score has to be 49.
# L,16,1. L=length, 16=constant, 1=multiplier of L (1xL). If you have a length of 33, then say 16+(33x1)=49.
# Matched bases=2 points. Mismatched bases=-6 points, and 2 base pair read gaps= -11 points
# (33x2)-6-11=49, therefore our alignment threshold allows 1 base pair mismatch and 1 gap, or 2 mismatches 

# Submit the job
qsub maps

# Create a list of the .sam files you made and check how many - should equal your amount of samples
ls -l *.sam | wc -l
ls *.sam > sams
cat sams | wc -l

# Next stage is compressing, sorting and indexing the SAM files, so they become BAM files:

# Create a file called s2b

>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

# Add module load samtools to top of s2b file

nano s2b

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N s2b # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o s2b.qlog # Name the file where to redirect standard output and error

module load samtools
samtools sort -O bam -o Pool_10_S10_L002_R1_001_ACCA.trim.bt2.bam Pool_10_S10_L002_R1_001_ACCA.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_ACCA.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_AGAC.trim.bt2.bam Pool_10_S10_L002_R1_001_AGAC.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_AGAC.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_AGTG.trim.bt2.bam Pool_10_S10_L002_R1_001_AGTG.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_AGTG.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_CATC.trim.bt2.bam Pool_10_S10_L002_R1_001_CATC.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_CATC.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_CTAC.trim.bt2.bam Pool_10_S10_L002_R1_001_CTAC.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_CTAC.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_GACT.trim.bt2.bam Pool_10_S10_L002_R1_001_GACT.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_GACT.trim.bt2.bam
samtools sort -O bam -o Pool_10_S10_L002_R1_001_GCTT.trim.bt2.bam Pool_10_S10_L002_R1_001_GCTT.trim.bt2.sam && samtools index Pool_10_S10_L002_R1_001_GCTT.trim.bt2.bam


# Run all commands listed in s2b file, add header as above
qsub s2b

# Do we have the correct number of bam files? Should be the same number as number of trim files
ls *.bam | wc -l  

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.


#===================================== ANGSD "Fuzzy" Genotyping ============================================
pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

# This first run-through of ANGSD is just to take a look at base qualities and coverage depth 
# We will run angsd again with filters informed by this first run-through

# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).

ls *.bam > bams
ls -l *.bam | wc -l
#Urgent: Before putting anything into R, you'll need a .csv file with the bam names and other data about each individual. If the order of bams is not the same as in the bams file on SCC, then samples will be mislabelled in R figures. It's a good idea to add any relevant metadata too

#------------------------------- Assessing base qualities and coverage depth -------------------------------------

# angsd settings:
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping =< 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -maxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd : the minimal number of individuals the site must be genotyped in. Reset to 50% of total N at this stage.

>angsdDD
nano angsdDD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsdDD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsdDD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -maxDepth 4660 -minInd 233"
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out dd 

qsub angsdDD


# Summarizing results (using modified script by Matteo Fumagalli)
module load R
Rscript /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/plotQC.R prefix=dd

# Proportion of sites covered at >5x:
cat quality.txt

# I've removed all low quality bams that have a proportion of sites covered at >5x under 0.4; the rest of the bams are in bams_all.csv
# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs 


#----------------------------------- Tech Rep Detection ------------------------------------------------------
# Looking at all samples including technical replicates to confirm accuracy of sequencing and also help determine filtering thresholds

## Detecting TRs (note: lower minind makes it easier to determine tech reps, but should be raised for subsequent analyses!)

>angsd_TR_all_mi50
nano angsd_TR_all_mi50

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_TR_all_mi50 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_TR_all_mi50.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 50 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams -GL 1 $FILTERS $TODO -P 1 -out angsd_TR_all_mi50

qsub angsd_TR_all_mi50

NSITES=`zcat angsd_TR_all_mi50.mafs.gz | wc -l` 
echo $NSITES
#31694 sites


#----------------------------------- Tech Rep Selection ------------------------------------------------------ 

cat quality.txt

# Now going to look at bams of tech reps to stick with the better of the two (better prop. of sites covered at >5x)

#268-G11 & 268-G22: Pool_10_Jaskiel_Lane2_S10_L004_R1_001_CTAC.trim.bt2.bam 0.691597161485455 > Pool_1_Jaskiel_Lane1_S1_L003_R1_001_AGAC.trim.bt2.bam 0.609519935159601
#281-G1: Pool_9_Jaskiel_Lane2_S9_L004_R1_001_TGGT.trim.bt2.bam 0.57977023819394 > Pool_11_Jaskiel_Lane2_S11_L004_R1_001_AGTG.trim.bt2.bam 0.391114605112956
#B11: Pool_8_Jaskiel_Lane1_S8_L003_R1_001_TGGT.trim.bt2.bam 0.931344870455589 > Pool_14_Jaskiel_Lane2_S14_L004_R1_001_TGGT.trim.bt2.bam 0.927000482594793 & Pool_15_S15_L002_R1_001_TCAC.trim.bt2.bam 0.900899577667227
#B12: Pool_14_Jaskiel_Lane2_S14_L004_R1_001_AGAC.trim.bt2.bam 0.827717637896629 & Pool_16_S16_L002_R1_001_TCAC.trim.bt2.bam 0.624755843810437 < Pool_8_Jaskiel_Lane1_S8_L003_R1_001_AGAC.trim.bt2.bam 0.848159228452337
#B13: Pool_15_Jaskiel_Lane2_S15_L004_R1_001_ACCA.trim.bt2.bam 0.854725311512307 < Pool_8_Jaskiel_Lane1_S8_L003_R1_001_ACCA.trim.bt2.bam 0.90211529814634
#B15: Pool_15_Jaskiel_Lane2_S15_L004_R1_001_CATC.trim.bt2.bam 0.813947662137126 < Pool_8_Jaskiel_Lane1_S8_L003_R1_001_CATC.trim.bt2.bam 0.877201538905538
#B16: Pool_8_Jaskiel_Lane1_S8_L003_R1_001_GTGA.trim.bt2.bam 0.914926127204862 > Pool_14_Jaskiel_Lane2_S14_L004_R1_001_GTGA.trim.bt2.bam 0.895720766491337
#B17: Pool_15_Jaskiel_Lane2_S15_L004_R1_001_TCAG.trim.bt2.bam 0.867398766004564 < Pool_8_Jaskiel_Lane1_S8_L003_R1_001_TCAG.trim.bt2.bam 0.905684095225631
#B23: Pool_15_Jaskiel_Lane2_S15_L004_R1_001_GACT.trim.bt2.bam 0.87548363026644 < Pool_8_Jaskiel_Lane1_S8_L003_R1_001_GACT.trim.bt2.bam 0.9170245151796
#Y20: Pool_15_Jaskiel_Lane2_S15_L004_R1_001_TGGT.trim.bt2.bam 0.946176516867479 & Pool_15_S15_L002_R1_001_GACT.trim.bt2.bam 0.856677165067073 < Pool_7_Jaskiel_Lane1_S7_L003_R1_001_TGGT.trim.bt2.bam 0.953680638244451
#Y21: Pool_15_Jaskiel_Lane2_S15_L004_R1_001_AGAC.trim.bt2.bam 0.80601536208809 < Pool_16_S16_L002_R1_001_GACT.trim.bt2.bam 0.8834173236762 > Pool_7_Jaskiel_Lane1_S7_L003_R1_001_AGAC.trim.bt2.bam 0.844064750630975
#Y25: Pool_15_Jaskiel_Lane2_S15_L004_R1_001_GTGA.trim.bt2.bam 0.911862499217695 < Pool_7_Jaskiel_Lane1_S7_L003_R1_001_GTGA.trim.bt2.bam 0.923540888153826
#268-G188: Pool_1_S1_L001_R1_001_TCAG.trim.bt2.bam 0.704454848705502 > Pool_8_S8_L001_R1_001_TCAG.trim.bt2.bam 0.622703106237599
#268-G66: Pool_2_S2_L001_R1_001_GCTT.trim.bt2.bam 0.721225379223316 > Pool_8_S8_L001_R1_001_GCTT.trim.bt2.bam 0.50270728203949
#268-G134: Pool_3_S3_L001_R1_001_CTAC.trim.bt2.bam 0.0941717912305791 & Pool_8_S8_L001_R1_001_CTAC.trim.bt2.bam 0.0853656795811435 <- GET RID OF BOTH
#268-G273: Pool_8_S8_L001_R1_001_TGTC.trim.bt2.bam 0.553758536784817 < Pool_4_S4_L001_R1_001_TGTC.trim.bt2.bam 0.567300116013541
#274-G50: Pool_5_S5_L001_R1_001_TCAC.trim.bt2.bam 0.584053755985297 > Pool_8_S8_L001_R1_001_TCAC.trim.bt2.bam 0.523210661096703
#281-G60: Pool_6_S6_L001_R1_001_GACT.trim.bt2.bam 0.898437010071253 > Pool_8_S8_L001_R1_001_GACT.trim.bt2.bam 0.895910531598728
#281-G203: Pool_8_S8_L001_R1_001_TGGT.trim.bt2.bam 0.750525547478687 > Pool_9_S9_L002_R1_001_TGGT.trim.bt2.bam 0.711435166858636
#281-G204: Pool_8_S8_L001_R1_001_AGAC.trim.bt2.bam 0.876131909932215 > Pool_9_S9_L002_R1_001_AGAC.trim.bt2.bam 0.683234122250635
#281-G205: Pool_8_S8_L001_R1_001_ACCA.trim.bt2.bam 0.872370903072389 > Pool_9_S9_L002_R1_001_ACCA.trim.bt2.bam 0.802791897826015
#281-G206: Pool_8_S8_L001_R1_001_AGTG.trim.bt2.bam 0.822907067849022 > Pool_9_S9_L002_R1_001_AGTG.trim.bt2.bam 0.598817170282787
#281-G209: Pool_8_S8_L001_R1_001_CATC.trim.bt2.bam 0.763366853793186 < Pool_9_S9_L002_R1_001_CATC.trim.bt2.bam 0.775322809868473
#281-G214: Pool_8_S8_L001_R1_001_GTGA.trim.bt2.bam 0.719138927986701 > Pool_9_S9_L002_R1_001_GTGA.trim.bt2.bam 0.688016199691803
#281-G246: Pool_11_S11_L002_R1_001_TGGT.trim.bt2.bam 0.676928858240945 > Pool_22_S6_L001_R1_001_TGGT.trim.bt2.bam 0.652200184339798
#281-G247: Pool_11_S11_L002_R1_001_AGAC.trim.bt2.bam 0.670736454464685 < Pool_22_S6_L001_R1_001_AGAC.trim.bt2.bam 0.717780184049497
#281-G248: Pool_11_S11_L002_R1_001_ACCA.trim.bt2.bam 0.862643696850512 > Pool_22_S6_L001_R1_001_ACCA.trim.bt2.bam 0.854173542985569
#281-G249: Pool_22_S6_L001_R1_001_AGTG.trim.bt2.bam 0.816192763204646 < Pool_11_S11_L002_R1_001_AGTG.trim.bt2.bam 0.85862891664411
#281-G250: Pool_11_S11_L002_R1_001_CATC.trim.bt2.bam 0.692983633333102 > Pool_22_S6_L001_R1_001_CATC.trim.bt2.bam 0.48088665951395
#281-G251: Pool_11_S11_L002_R1_001_GTGA.trim.bt2.bam 0.757257364750595 > Pool_22_S6_L001_R1_001_GTGA.trim.bt2.bam 0.688183311636949
#281-G233: Pool_10_S10_L002_R1_001_ACCA.trim.bt2.bam 0.647753676571331 < Pool_16_S16_L002_R1_001_ACCA.trim.bt2.bam 0.656886239812877
#281-G237: Pool_10_S10_L002_R1_001_CATC.trim.bt2.bam 0.894635823913688 > Pool_16_S16_L002_R1_001_CATC.trim.bt2.bam 0.88273364801248
#281-G252: Pool_11_S11_L002_R1_001_TCAG.trim.bt2.bam 0.693736108815314 < Pool_16_S16_L002_R1_001_TCAG.trim.bt2.bam 0.780715479631195
#268-G194: Pool_13_S13_L002_R1_001_CTAC.trim.bt2.bam 0.627694653765432 > Pool_16_S16_L002_R1_001_CTAC.trim.bt2.bam 0.454416585433252
#268-G248: Pool_14_S14_L002_R1_001_TGTC.trim.bt2.bam 0.825628377211454 > Pool_16_S16_L002_R1_001_TGTC.trim.bt2.bam 0.765302837931228
#P7: Pool_22_S6_L001_R1_001_TCAG.trim.bt2.bam 0.122572714427212 > Pool_23_S7_L001_R1_001_TCAG.trim.bt2.bam 0.118404084287964 <-- GET RID OF BOTH?
#P8: Pool_22_S6_L001_R1_001_GCTT.trim.bt2.bam 0.51087003160051 > Pool_23_S7_L001_R1_001_GCTT.trim.bt2.bam 0.499552676343243
#P9: Pool_22_S6_L001_R1_001_CTAC.trim.bt2.bam 0.876802860530831 < Pool_23_S7_L001_R1_001_CTAC.trim.bt2.bam 0.884488588314413
#P10: Pool_22_S6_L001_R1_001_TGTC.trim.bt2.bam 0.83930787517588 < Pool_23_S7_L001_R1_001_TGTC.trim.bt2.bam 0.850931171573292
#P11: Pool_22_S6_L001_R1_001_TCAC.trim.bt2.bam 0.851163776437385 < Pool_23_S7_L001_R1_001_TCAC.trim.bt2.bam 0.858114408882457
#P12: Pool_22_S6_L001_R1_001_GACT.trim.bt2.bam 0.839345068088488 < Pool_23_S7_L001_R1_001_GACT.trim.bt2.bam 0.851070545749508
#261-G21: Pool_10_Jaskiel_Lane2_S10_L004_R1_001_AGAC.trim.bt2.bam 0.70846361862503 > Pool_15_S15_L002_R1_001_TGGT.trim.bt2.bam 0.595432006996303
#261-G22: Pool_10_Jaskiel_Lane2_S10_L004_R1_001_ACCA.trim.bt2.bam 0.838460423436809 > Pool_16_S16_L002_R1_001_TGGT.trim.bt2.bam 0.709792651097786
#261-G23: Pool_4_Jaskiel_Lane1_S4_L003_R1_001_TGGT.trim.bt2.bam 0.942801048063264 > Pool_15_S15_L002_R1_001_AGAC.trim.bt2.bam 0.894785835500684
#268-G6: Pool_9_Jaskiel_Lane2_S9_L004_R1_001_TCAG.trim.bt2.bam 0.71353266553661 < Pool_15_S15_L002_R1_001_ACCA.trim.bt2.bam 0.727858500775601
#268-G18: Pool_15_S15_L002_R1_001_AGTG.trim.bt2.bam 0.745570184650985 < Pool_2_Jaskiel_Lane1_S2_L003_R1_001_TGGT.trim.bt2.bam 0.855701713386485
#274-G6: Pool_10_Jaskiel_Lane2_S10_L004_R1_001_TGGT.trim.bt2.bam 0.73431748641204 > Pool_15_S15_L002_R1_001_CATC.trim.bt2.bam 0.649633236724969
#274-G12: Pool_1_Jaskiel_Lane1_S1_L003_R1_001_AGTG.trim.bt2.bam 0.493800510813582 < Pool_15_S15_L002_R1_001_GTGA.trim.bt2.bam 0.753822781352698
#281-G6: Pool_15_S15_L002_R1_001_TCAG.trim.bt2.bam 0.252977242495217 < Pool_11_Jaskiel_Lane2_S11_L004_R1_001_CTAC.trim.bt2.bam 0.47914619092284
#281-G22: Pool_15_S15_L002_R1_001_GCTT.trim.bt2.bam 0.522390769026508 < Pool_4_Jaskiel_Lane1_S4_L003_R1_001_CTAC.trim.bt2.bam 0.523332932692992
#287-G8: Pool_10_Jaskiel_Lane2_S10_L004_R1_001_GCTT.trim.bt2.bam 0.548975954884898 > Pool_15_S15_L002_R1_001_CTAC.trim.bt2.bam 0.463369865699923
#287-G13: Pool_15_S15_L002_R1_001_TGTC.trim.bt2.bam 0.422215126537785 < Pool_9_Jaskiel_Lane2_S9_L004_R1_001_TCAC.trim.bt2.bam 0.454008351855823
#261-G26: Pool_10_Jaskiel_Lane2_S10_L004_R1_001_CATC.trim.bt2.bam 0.779633978348889 < Pool_16_S16_L002_R1_001_AGAC.trim.bt2.bam 0.783776966104911
#268-G21: Pool_16_S16_L002_R1_001_AGTG.trim.bt2.bam 0.700968726402163 < Pool_1_Jaskiel_Lane1_S1_L003_R1_001_TGGT.trim.bt2.bam 0.729453258440989
#274-G14: Pool_4_Jaskiel_Lane1_S4_L003_R1_001_GTGA.trim.bt2.bam 0.807561496865569 > Pool_16_S16_L002_R1_001_GTGA.trim.bt2.bam 0.747946642944272
#281-G23: Pool_16_S16_L002_R1_001_GCTT.trim.bt2.bam 0.694371029455565 > Pool_2_Jaskiel_Lane1_S2_L003_R1_001_AGAC.trim.bt2.bam 0.612409304115472
#268-G17 & 268-G266: Pool_14_Jaskiel_Lane2_S14_L004_R1_001_TCAC.trim.bt2.bam 0.269579754496615 < Pool_4_S4_L001_R1_001_CTAC.trim.bt2.bam 0.867586087916686
#268-G14 & 268-G251: Pool_14_Jaskiel_Lane2_S14_L004_R1_001_GCTT.trim.bt2.bam 0.278375903853978 < Pool_14_S14_L002_R1_001_TCAC.trim.bt2.bam 0.881108758519965
#268-G13 & 268-G18: Pool_14_Jaskiel_Lane2_S14_L004_R1_001_TGTC.trim.bt2.bam 0.613516556461277 < Pool_2_Jaskiel_Lane1_S2_L003_R1_001_TGGT.trim.bt2.bam 0.855701713386485
#268-G16 & 268-G261: Pool_13_Jaskiel_Lane2_S13_L004_R1_001_TGGT.trim.bt2.bam 0.541977248171747 < Pool_4_S4_L001_R1_001_GCTT.trim.bt2.bam 0.830989439241658

# I've removed all low quality bams that are under 0.4 as well as TRs; the rest of the bams are in bams_all.csv

#Also removing 3 weird "noisy" samples, 268-G273, 281-G48, and 287-G53
#Pool_4_S4_L001_R1_001_TGTC.trim.bt2.bam, Pool_2_Jaskiel_Lane1_S2_L003_R1_001_GTGA.trim.bt2.bam, Pool_12_S12_L002_R1_001_TGTC.trim.bt2.bam


# - - - - - - - - - - - - - - - - - All samples, no TRs or low coverage samples - - - - - - - - - - - - - - - - - -

#Filtering out sites in linkage (can drive population structure)

# filtering sites to work on - use only filters that do not distort allele frequency
# set minInd to 75-90% of the total number fo individuals in the project
# if you are doing any other RAD than 2bRAD or GBS, remove '-sb_pval 1e-5' from FILTERS
>angsd_AllSites
nano angsd_AllSites

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_AllSites # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_AllSites.qlog # Name the file where to redirect standard output and error

module load angsd

FILTERS="-uniqueOnly 1 -remove_bads 1  -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 289"
TODO="-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2 -doGeno 8"

angsd -b bams_all_noTR -GL 1 $FILTERS $TODO -P 1 -out AllSites

qsub angsd_AllSites

NSITES=`zcat AllSites.mafs.gz | wc -l` 
echo $NSITES
#704405

# Collecting and indexing filter-passing sites

zcat AllSites.mafs.gz | cut -f 1,2 | tail -n +2 >AllSites

>indexing_AllSites
nano indexing_AllSites

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N indexing_AllSites # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o indexing_AllSites.qlog # Name the file where to redirect standard output and error

module load angsd
angsd sites index AllSites

qsub indexing_AllSites

zcat AllSites.mafs.gz | tail -n +2 | cut -f 1,2 > mc1.sites

##Removing linked sites##
#Run as a job

>ngsld_unlinked_5kb
nano ngsld_unlinked_5kb

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsld_unlinked_5kb # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsld_unlinked_5kb.qlog # Name the file where to redirect standard output and error

module load ngsld
NS=`zcat AllSites.geno.gz| wc -l`
NB=`cat bams_all_noTR | wc -l`
ngsLD --geno AllSites.geno.gz --probs 1 --n_ind $NB --n_sites $NS --max_kb_dist 5 --pos AllSites --out AllSites_5kb.LD --n_threads 1 --extend_out 1 --min_maf 0.05

qsub ngsld_unlinked_5kb

>prune_LD_5kb
nano prune_LD_5kb

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N prune_LD_5kb # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o prune_LD_5kb.qlog # Name the file where to redirect standard output and error

module load perl
module load ngsld
cat AllSites_5kb.LD | cut -f 1,3,5- | perl /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/prune_graph.pl --max_kb_dist 5 --min_weight 0.2 --weight_type a > 0.2_unlinked

qsub prune_LD_5kb

sed 's/:/\t/g' 0.2_unlinked > 0.2_unlinked.sites
awk  '$2!=""' 0.2_unlinked.sites > 0.2_unlinked.sites.tmp; mv 0.2_unlinked.sites.tmp 0.2_unlinked.sites
sort -k1 0.2_unlinked.sites > 0.2_unlinked.sites_2.tmp; mv 0.2_unlinked.sites_2.tmp 0.2_unlinked.sites

angsd sites index 0.2_unlinked.sites

>final_sites_5kb
nano final_sites_5kb

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N final_sites_5kb # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o final_sites_5kb.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -hetbias_pval 1e-5 -minInd 289 -sb_pval 1e-5"
TODO='-doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 2'
angsd -sites 0.2_unlinked.sites -b bams_all_noTR -GL 1 $FILTERS $TODO -P 1 -out 0.2_unlinked.sites

qsub final_sites_5kb

NSITES=`zcat 0.2_unlinked.sites.mafs.gz | wc -l`
echo $NSITES
#197052

zcat 0.2_unlinked.sites.mafs.gz | cut -f 1,2 | tail -n +2 > finalsites_0.2_unlinked
angsd sites index finalsites_0.2_unlinked


#Without LD pruning 
# minind 289 (80% of 361) minQ 25 w/ -setMinDepthInd 5

>angsd_all_mid5
nano angsd_all_mid5

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_all_mid5 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_all_mid5.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 289 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_all_noTR -GL 1 $FILTERS $TODO -P 1 -out angsd_all_mid5

qsub angsd_all_mid5

NSITES=`zcat angsd_all_mid5.mafs.gz | wc -l` 
echo $NSITES
#2777 sites


#Without linked sites
>angsd_all_mid5_noLD
nano angsd_all_mid5_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_all_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_all_mid5_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 289 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_all_noTR -sites finalsites_0.2_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_all_mid5_noLD

qsub angsd_all_mid5_noLD

NSITES=`zcat angsd_all_mid5_noLD.mafs.gz | wc -l` 
echo $NSITES
#595 sites


>angsd_all_noLD
nano angsd_all_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_all_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_all_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 289 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_all_noTR -sites finalsites_0.2_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_all_noLD

qsub angsd_all_noLD

NSITES=`zcat angsd_all_noLD.mafs.gz | wc -l` 
echo $NSITES
#3200 sites


#No Auxis (357 so 80% is 286)
>angsd_all_mid5_noLD_nofrig
nano angsd_all_mid5_noLD_nofrig

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_all_mid5_noLD_nofrig # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_all_mid5_noLD_nofrig.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 286 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_all_noauxis -sites finalsites_0.2_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_all_mid5_noLD_nofrig

qsub angsd_all_mid5_noLD_nofrig

NSITES=`zcat angsd_all_mid5_noLD_nofrig.mafs.gz | wc -l` 
echo $NSITES
#610

#================================= NGSadmix and ADMIXTURE ========================================
pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

# First, install:
wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

# --- All taxa no technical replicates or genotyping duplicates (normal filters)
# NgsAdmix for K from 1 to 10 : do not run if the dataset contains clones or genotyping replicates!

>ngsadmix_all_mid5_noLD
nano ngsadmix_all_mid5_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_all_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_all_mid5_noLD.qlog # Name the file where to redirect standard output and error

for K in `seq 1 10` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_all_mid5_noLD.beagle.gz -K $K -P 1 -o ngsadmix_all_mid5_noLD_k${K};
done

qsub ngsadmix_all_mid5_noLD

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot

##Trying ADMIXTURE to find optimal K for all taxa, but not sure if it'll work because it's using called SNPs
#BCF -> VCF
module load htslib/1.16
module load bcftools/1.16

bcftools convert -O v -o angsd_all_mid5_noLD.vcf angsd_all_mid5_noLD.bcf

module load admixture/1.3.0
module load plink/1.90b6.4

#all_mid5_noLD
plink --vcf angsd_all_mid5_noLD.vcf --make-bed --allow-extra-chr 0 --out angsd_all_mid5_noLD --const-fid 0
for K in `seq 1 10`; \
do admixture --cv angsd_all_mid5_noLD.bed $K | tee myresult_all_mid5_noLD_${K}.out; done

# use this to check K of least CV error:
grep -h CV myresult_all_mid5_noLD_*.out

CV error (K=10): 0.22820
CV error (K=1): 0.90770
CV error (K=2): 0.25892
CV error (K=3): 0.26050
CV error (K=4): 0.22362
CV error (K=5): 0.22691
CV error (K=6): 0.21598 <- lowest
CV error (K=7): 0.21963
CV error (K=8): 0.21896
CV error (K=9): 0.22520

#Note: to get the least CV error, you either need to rename chr. (hard) or use --allow-extra-chr 0 
#So run with --allow-extra-chr 0 first then run again without for PCAdapt

#Without Auxis
>ngsadmix_all_mid5_noLD_nofrig
nano ngsadmix_all_mid5_noLD_nofrig

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_all_mid5_noLD_nofrig # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_all_mid5_noLD_nofrig.qlog # Name the file where to redirect standard output and error

for K in `seq 1 10` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_all_mid5_noLD_nofrig.beagle.gz -K $K -P 1 -o ngsadmix_all_mid5_noLD_nofrig_k${K};
done

qsub ngsadmix_all_mid5_noLD_nofrig

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot

##Trying ADMIXTURE to find optimal K for all taxa, but not sure if it'll work because it's using called SNPs
#BCF -> VCF
module load htslib/1.16
module load bcftools/1.16

bcftools convert -O v -o angsd_all_mid5_noLD_nofrig.vcf angsd_all_mid5_noLD_nofrig.bcf

module load admixture/1.3.0
module load plink/1.90b6.4

#all_mid5_noLD_nofrig
plink --vcf angsd_all_mid5_noLD_nofrig.vcf --make-bed --allow-extra-chr 0 --out angsd_all_mid5_noLD_nofrig --const-fid 0
for K in `seq 1 10`; \
do admixture --cv angsd_all_mid5_noLD_nofrig.bed $K | tee myresult_all_mid5_noLD_nofrig_${K}.out; done

# use this to check K of least CV error:
grep -h CV myresult_all_mid5_noLD_nofrig_*.out
CV error (K=10): 0.22639
CV error (K=1): 0.91182
CV error (K=2): 0.25597
CV error (K=3): 0.25978
CV error (K=4): 0.22295
CV error (K=5): 0.22561
CV error (K=6): 0.21536 <- lowest still 
CV error (K=7): 0.21818
CV error (K=8): 0.22049
CV error (K=9): 0.22230


#Skipjack (ref)
>ngsadmix_skj_mid5_noLD
nano ngsadmix_skj_mid5_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_skj_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_skj_mid5_noLD.qlog # Name the file where to redirect standard output and error

for K in `seq 1 10` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_skj_mid5_noLD.beagle.gz -K $K -P 1 -o ngsadmix_skj_mid5_noLD_k${K};
done

qsub ngsadmix_skj_mid5_noLD

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot

##Trying ADMIXTURE to find optimal K for all taxa
#BCF -> VCF
module load htslib/1.16
module load bcftools/1.16

bcftools convert -O v -o angsd_skj_mid5_noLD.vcf angsd_skj_mid5_noLD.bcf

module load admixture/1.3.0
module load plink/1.90b6.4

plink --vcf angsd_skj_mid5_noLD.vcf --make-bed --allow-extra-chr 0 --out angsd_skj_mid5_noLD --const-fid 0
for K in `seq 1 10`; \
do admixture --cv angsd_skj_mid5_noLD.bed $K | tee myresult_skj_mid5_noLD_${K}.out; done

# use this to check K of least CV error:
grep -h CV myresult_skj_mid5_noLD_*.out

CV error (K=10): 0.58352
CV error (K=1): 0.50759 <- lowest
CV error (K=2): 0.52303
CV error (K=3): 0.53507
CV error (K=4): 0.54495
CV error (K=5): 0.54791
CV error (K=6): 0.55230
CV error (K=7): 0.56591
CV error (K=8): 0.57169
CV error (K=9): 0.57592


#Yellowfin
>ngsadmix_yft_mid5_noLD
nano ngsadmix_yft_mid5_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_yft_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_yft_mid5_noLD.qlog # Name the file where to redirect standard output and error

for K in `seq 1 10` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_yft_mid5_noLD.beagle.gz -K $K -P 1 -o ngsadmix_yft_mid5_noLD_k${K};
done

qsub ngsadmix_yft_mid5_noLD

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot

##Trying ADMIXTURE to find optimal K for all taxa
#BCF -> VCF
module load htslib/1.16
module load bcftools/1.16

bcftools convert -O v -o angsd_yft_mid5_noLD.vcf angsd_yft_mid5_noLD.bcf

module load admixture/1.3.0
module load plink/1.90b6.4

plink --vcf angsd_yft_mid5_noLD.vcf --make-bed --allow-extra-chr 0 --out angsd_yft_mid5_noLD --const-fid 0
for K in `seq 1 10`; \
do admixture --cv angsd_yft_mid5_noLD.bed $K | tee myresult_yft_mid5_noLD_${K}.out; done

# use this to check K of least CV error:
grep -h CV myresult_yft_mid5_noLD_*.out

CV error (K=10): 0.58150
CV error (K=1): 0.61330
CV error (K=2): 0.40382 <- lowest
CV error (K=3): 0.42535
CV error (K=4): 0.44700
CV error (K=5): 0.46651
CV error (K=6): 0.47996
CV error (K=7): 0.50965
CV error (K=8): 0.52031
CV error (K=9): 0.54739


#Bigeye
>ngsadmix_bet_mid5_noLD
nano ngsadmix_bet_mid5_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix_bet_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o ngsadmix_bet_mid5_noLD.qlog # Name the file where to redirect standard output and error

for K in `seq 1 10` ;
do
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/NGSadmix -likes angsd_bet_mid5_noLD.beagle.gz -K $K -P 1 -o ngsadmix_bet_mid5_noLD_k${K};
done

qsub ngsadmix_bet_mid5_noLD

# scp *qopt files to local machine, use admixturePlotting_v5.R to plot

##Trying ADMIXTURE to find optimal K for all taxa
#BCF -> VCF
module load htslib/1.16
module load bcftools/1.16

bcftools convert -O v -o angsd_bet_mid5_noLD.vcf angsd_bet_mid5_noLD.bcf

module load admixture/1.3.0
module load plink/1.90b6.4

plink --vcf angsd_bet_mid5_noLD.vcf --make-bed --allow-extra-chr 0 --out angsd_bet_mid5_noLD --const-fid 0
for K in `seq 1 10`; \
do admixture --cv angsd_bet_mid5_noLD.bed $K | tee myresult_bet_mid5_noLD_${K}.out; done

# use this to check K of least CV error:
grep -h CV myresult_bet_mid5_noLD_*.out

CV error (K=10): 1.08708
CV error (K=1): 0.52436 <- lowest
CV error (K=2): 0.57006
CV error (K=3): 0.64617
CV error (K=4): 0.68057
CV error (K=5): 0.75060
CV error (K=6): 0.78820
CV error (K=7): 0.86285
CV error (K=8): 0.90836
CV error (K=9): 1.00806

#================================= Demographic Analyses ========================================

## SFS work with individual species groups

#Skipjack (203 individuals, so 80% is 162)
>sfs_skj_dem
nano sfs_skj_dem

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_dem # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_dem.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_skj -GL 1 -P 1 -minInd 162 $TODO -out skj_dem

qsub sfs_skj_dem


#Yellowfin 
>sfs_yft_dem
nano sfs_yft_dem

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_dem # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_dem.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft -GL 1 -P 1 -minInd 77 $TODO -out yft_dem

qsub sfs_yft_dem


#Bigeye
>sfs_bet_dem
nano sfs_bet_dem

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_dem # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_dem.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_bet -GL 1 -P 1 -minInd 46 $TODO -out bet_dem

qsub sfs_bet_dem


# generating per-population SFS
realSFS skj_dem.saf.idx >skj_dem.sfs
realSFS yft_dem.saf.idx >yft_dem.sfs
realSFS bet_dem.saf.idx >bet_dem.sfs

# writing down 2d-SFS priors - SKJ vs YFT 
realSFS skj_dem.saf.idx yft_dem.saf.idx -P 24 > skj_yft_dem.sfs ; realSFS fst index skj_dem.saf.idx yft_dem.saf.idx -sfs skj_yft_dem.sfs -fstout skj_yft_dem

# global Fst between populations (after LD pruning)
realSFS fst stats skj_yft_dem.fst.idx

#output:
FST.Unweight[nObs:143797]:0.018354 Fst.Weight:0.767844

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# writing down 2d-SFS priors - SKJ vs BET 
realSFS skj_dem.saf.idx bet_dem.saf.idx -P 24 > skj_bet_dem.sfs ; realSFS fst index skj_dem.saf.idx bet_dem.saf.idx -sfs skj_bet_dem.sfs -fstout skj_bet_dem

# global Fst between populations (after LD pruning)
realSFS fst stats skj_bet_dem.fst.idx

#output:
FST.Unweight[nObs:152273]:0.020296 Fst.Weight:0.786622


# writing down 2d-SFS priors - YFT vs BET 
realSFS yft_dem.saf.idx bet_dem.saf.idx -P 24 > yft_bet_dem.sfs ; realSFS fst index yft_dem.saf.idx bet_dem.saf.idx -sfs yft_bet_dem.sfs -fstout yft_bet_dem

# global Fst between populations (after LD pruning)
realSFS fst stats yft_bet_dem.fst.idx

#output:
FST.Unweight[nObs:145287]:0.007210 Fst.Weight:0.305582


#-----------------------------------------Obtaining bootstrapped SFS with ANGSD-------------------------------------------
## Analysis in moments (see: https://github.com/z0on/AFS-analysis-with-moments/blob/master/README.md)
## Pairwise bootstrapping SFS

#Now we generate the bootstrapped data (100 series of 5 bootstraps) for each pairwise comparison
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta

#starting off with SKJ and YFT
>skj_yft_b100
for B in `seq 1 100`; do
echo "realSFS skj_dem.saf.idx yft_dem.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 5 -P 1 -resample_chr 1 >skj_yft_boot_$B">>skj_yft_b100;
done

>skj_yft_b100_array
nano skj_yft_b100_array

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N skj_yft_b100_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m a
#$ -pe omp 1
#$ -t 1-100
#$ -j y # Join standard output and error to a single file
#$ -o skj_yft_b100_array.qlog # Name the file where to redirect standard output and error

module load angsd
readarray -t CMD_ARRAY < skj_yft_b100
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub skj_yft_b100_array


#SKJ and BET
>skj_bet_b100
for B in `seq 1 100`; do
echo "realSFS skj_dem.saf.idx bet_dem.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 5 -P 1 -resample_chr 1 >skj_bet_boot_$B">>skj_bet_b100;
done

>skj_bet_b100_array
nano skj_bet_b100_array

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N skj_bet_b100_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m a
#$ -pe omp 1
#$ -t 1-100
#$ -j y # Join standard output and error to a single file
#$ -o skj_bet_b100_array.qlog # Name the file where to redirect standard output and error

module load angsd
readarray -t CMD_ARRAY < skj_bet_b100
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub skj_bet_b100_array


#YFT and BET
>yft_bet_b100
for B in `seq 1 100`; do
echo "realSFS yft_dem.saf.idx bet_dem.saf.idx -ref $GENOME_REF -anc $GENOME_REF -bootstrap 5 -P 1 -resample_chr 1 >yft_bet_boot_$B">>yft_bet_b100;
done

>yft_bet_b100_array
nano yft_bet_b100_array

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N yft_bet_b100_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m a
#$ -pe omp 1
#$ -t 1-100
#$ -j y # Join standard output and error to a single file
#$ -o yft_bet_b100_array.qlog # Name the file where to redirect standard output and error

module load angsd
readarray -t CMD_ARRAY < yft_bet_b100
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub yft_bet_b100_array


#Computing SFS dimensions
#SKJ
export N1=`wc -l bams_skj | cut -f 1 -d " "`
export NG1=`echo "($N1*2)+1" | bc`

#YFT
export N2=`wc -l bams_yft | cut -f 1 -d " "`
export NG2=`echo "($N2*2)+1" | bc`

#BET
export N3=`wc -l bams_bet | cut -f 1 -d " "`
export NG3=`echo "($N3*2)+1" | bc`

#Averaging 5-bootstrap batches
#SKJ and YFT
for B in `seq 1 100`; do
echo "$NG1 $NG2" >skj_yft_boot_${B}.sfs;
cat skj_yft_boot_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> skj_yft_boot_${B}.sfs;
done

#Move to its own folder
mv skj_yft_boot_* skj_yft_moments_analysis/


#SKJ and BET
for B in `seq 1 100`; do
echo "$NG1 $NG3" >skj_bet_boot_${B}.sfs;
cat skj_bet_boot_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> skj_bet_boot_${B}.sfs;
done

#Move to its own folder
mv skj_bet_boot_* skj_bet_moments_analysis/


#YFT and BET
for B in `seq 1 100`; do
echo "$NG2 $NG3" >yft_bet_boot_${B}.sfs;
cat yft_bet_boot_${B} | awk '{for (i=1;i<=NF;i++){a[i]+=$i;}} END {for (i=1;i<=NF;i++){printf "%.3f", a[i]/NR; printf "\t"};printf "\n"}' >> yft_bet_boot_${B}.sfs;
done

#Move to its own folder
mv yft_bet_boot_* yft_bet_moments_analysis/


#FLAG: might end up excluding moments, in which case delete this. If not, finish running it after figuring out SU issue, then update R doc to have the files used for this analysis (all others done)
#-------------------------------- AFS Analysis in Moments ---------------------------------
#Continues on the above code from https://github.com/z0on/AFS-analysis-with-moments

#First of all, install moments. The example below would clone it into the user's home directory and install it for a specific user.
cd
pwd
/usr3/graduate/jaskielj

git clone --branch devel https://bitbucket.org/simongravel/moments.git 
cd moments
module load python3
python3 setup.py build_ext --inplace

# add moments to $PYTHONPATH (consider adding this line to your .bashrc):
export PYTHONPATH=$PYTHONPATH:/usr3/graduate/jaskielj/moments

cd 
git clone https://github.com/z0on/AFS-analysis-with-moments.git
cd AFS-analysis-with-moments
mkdir work
cp multimodel_inference/py3_v1/* work/

#NOTE: all code examples here assume the repository is cloned in the home directory, ~/. If you cloned it elsewhere, make sure to replace ~/AFS-analysis-with-moments in all examples with the actual path.
#Overview:
#The first step is model selection, where we run all possible models on 10 bootstrapped SFS. We run each model on each bootstrap six times (six random restarts), to make sure the model converges to its best likelihood at least once. All these commands are written by the R script modSel_write.R. Then we use the R script modSel_summary.R to select the best-fitted instance (out of 6) for each model for each bootstrap, and compare the AIC scores for all models. The best model is the one with the lowest median AIC score among bootstrap replicates.
#The second step is running the winning model on 100 bootstrapped SFS, to evaluate parameter uncertainties. The commands for this stage are actually written by the modSel_summary.R script. Once again, we are doing 6 random restarts for each bootstrap replicate. The parameter meanings and uncertainties are deciphered by the third R script that we have, bestBoot_summary.R. All three R scripts are designed for command-line usage.
#NOTE for version 1: the model names are not fully standardized to the above convention, please see work/moments_multimodels.xlsx for their structure.

#Basic example:
Rscript ~/AFS-analysis-with-moments/work/modSel_write.R contrast=p12 args="p1 p2 16 16 0.02 0.005"
#Where...
#contrast : the name of population comparison. It should match the leading part of the bootstapped SFS names (in example here, p12)
#args : list of parameters for AFS models, in the following order: 
#name of pop1, name of pop2 
#projection for pop1, projection for pop2
#mutation rate (per genotyped portion of the genome per generation)
#generation time in thousands of years

#SKJ: 203 x 1.6 = 325
#YFT: 96 x 1.6 = 154
#BET: 57 x 1.6 = 91

#mutation rate per base should be multiplied by the number of sites in the portion of the genome under analysis 
#Where L is the number of monomorphic and polymorphic sites used for the moments analysis (can be extracted from the number of lines in your .maf file) 
#Mutation rates and generation times:
#YFT: mutation rate = 7.3x10^-9 per site per generation (Barth et al 2017)

#Generation times:
#SKJ: 1 yr
#YFT: 2 yrs
#BET: 3 yrs
#Avg: 2 yrs

#sites used in moments analysis = 197052 SNPs
#mu = 7.3x10^-9
#so mutations per genotyped portion of the genome per generation = 7.3x10^-9 x 197052 = 0.0014


#Model Selection - SKJ and YFT
cd /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/skj_yft_moments_analysis

>bestmod_skj_yft
nano bestmod_skj_yft

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bestmod_skj_yft # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o bestmod_skj_yft.qlog # Name the file where to redirect standard output and error

module load R
Rscript ~/AFS-analysis-with-moments/work/modSel_write.R contrast=skj_yft_boot args="skj yft 325 154 0.0014 0.002"

qsub bestmod_skj_yft

#Commands get printed to skj_yft_boot.modsel.runs

#Now that we have our file with commands in it, we'll make another file that has the task array code that references the command file
>skj_yft_moments_array
nano skj_yft_moments_array

#add the header and load python, add PYTHONPATH
#remove the #$ -m bea option, or you will have a very bad time

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N skj_yft_moments_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -pe omp 1
#$ -t 1-6480
#$ -j y # Join standard output and error to a single file
#$ -o skj_yft_moments_array.qlog # Name the file where to redirect standard output and error

module load python3
export PYTHONPATH=$PYTHONPATH:/usr3/graduate/jaskielj/moments
readarray -t CMD_ARRAY < skj_yft_boot.modsel.runs
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub skj_yft_moments_array

#so everything get printed to skj_yft_boot.modsel

#See modSel_summary_JJ.R for my custom script based on James' modification of Misha's script
#FLAG: Add to main R code file

sed -i 's/\x00//g' skj_yft_boot.modsel
grep -a RESULT skj_yft_boot.modsel -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT  >skj_yft_boot.modsel.res1
#go through .res1 and make sure there aren't any lines that were cut off (delete those lines)
#Go through modSel_summary_skj_yft.R line by line; it doesn't really work if you run it as a job

NREPS=6 # number of random restarts per model per bootstrap rep
>winner.mods
for i in `seq 1 $NREPS`;do
awk '{print $2}' *.winmod | sed 's/$/.py/'  >>winner.mods;
done

sed -i -e 's/^/\/usr3\/graduate\/jaskielj\/AFS-analysis-with-moments\/work\//' winner.mods

CONTRAST=skj_yft_boot
ARGS="skj yft 325 154 0.0014 0.002"

>winner.modsel
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat winner.mods | wc -l`
>winner.args
>${CONTRAST}.winboots
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS >>${CONTRAST}.winboots" >>winner.args;
done;
paste winner.mods winner.args -d " " >>winner.modsel;
done

#Running the winning model as an array
>skj_yft_winmodsel_array
nano skj_yft_winmodsel_array

#add the header and load python, add PYTHONPATH

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N skj_yft_winmodsel_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -pe omp 1
#$ -t 1-600
#$ -j y # Join standard output and error to a single file
#$ -o skj_yft_winmodsel_array.qlog # Name the file where to redirect standard output and error

module load python3
export PYTHONPATH=$PYTHONPATH:/usr3/graduate/jaskielj/moments
readarray -t CMD_ARRAY < winner.modsel
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub skj_yft_winmodsel_array


#Model Selection - SKJ and BET
cd /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/skj_bet_moments_analysis

>bestmod_skj_bet
nano bestmod_skj_bet

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bestmod_skj_bet # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o bestmod_skj_bet.qlog # Name the file where to redirect standard output and error

module load R
Rscript ~/AFS-analysis-with-moments/work/modSel_write.R contrast=skj_bet_boot args="skj bet 325 91 0.0014 0.002"

qsub bestmod_skj_bet

#Commands get printed to skj_bet_boot.modsel.runs
nano skj_bet_boot.modsel.runs

#Now that we have our file with commands in it, we'll make another file that has the task array code that references the command file
>skj_bet_moments_array
nano skj_bet_moments_array

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N skj_bet_moments_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -pe omp 1
#$ -t 1-6480
#$ -j y # Join standard output and error to a single file
#$ -o skj_bet_moments_array.qlog # Name the file where to redirect standard output and error

module load python3
export PYTHONPATH=$PYTHONPATH:/usr3/graduate/jaskielj/moments
readarray -t CMD_ARRAY < skj_bet_boot.modsel.runs
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub skj_bet_moments_array

#so everything get printed to skj_bet_boot.modsel

sed -i 's/\x00//g' skj_bet_boot.modsel
grep -a RESULT skj_bet_boot.modsel -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT  >skj_bet_boot.modsel.res1
#go through .res1 and make sure there aren't any lines that were cut off (delete those lines)
#Go through modSel_summary_skj_bet.R line by line; it doesn't really work if you run it as a job

NREPS=6 # number of random restarts per model per bootstrap rep
>winner.mods
for i in `seq 1 $NREPS`;do
awk '{print $2}' *.winmod | sed 's/$/.py/'  >>winner.mods;
done

sed -i -e 's/^/\/usr3\/graduate\/jaskielj\/AFS-analysis-with-moments\/work\//' winner.mods

CONTRAST=skj_bet_boot
ARGS="skj bet 325 91 0.0014 0.002"

>winner.modsel
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat winner.mods | wc -l`
>winner.args
>${CONTRAST}.winboots
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS >>${CONTRAST}.winboots" >>winner.args;
done;
paste winner.mods winner.args -d " " >>winner.modsel;
done

#Running the winning model as an array
>skj_bet_winmodsel_array
nano skj_bet_winmodsel_array

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N skj_bet_winmodsel_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -pe omp 1
#$ -t 1-600
#$ -j y # Join standard output and error to a single file
#$ -o skj_bet_winmodsel_array.qlog # Name the file where to redirect standard output and error

module load python3
export PYTHONPATH=$PYTHONPATH:/usr3/graduate/jaskielj/moments
readarray -t CMD_ARRAY < winner.modsel
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub skj_bet_winmodsel_array


#Model Selection - YFT and BET
cd /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/yft_bet_moments_analysis

>bestmod_yft_bet
nano bestmod_yft_bet

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bestmod_yft_bet # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o bestmod_yft_bet.qlog # Name the file where to redirect standard output and error

module load R
Rscript ~/AFS-analysis-with-moments/work/modSel_write.R contrast=yft_bet_boot args="yft bet 154 91 0.0014 0.002"

qsub bestmod_yft_bet

#Commands get printed to yft_bet_boot.modsel.runs

#Now that we have our file with commands in it, we'll make another file that has the task array code that references the command file
>yft_bet_moments_array
nano yft_bet_moments_array

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N yft_bet_moments_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -pe omp 1
#$ -t 1-6480
#$ -j y # Join standard output and error to a single file
#$ -o yft_bet_moments_array.qlog # Name the file where to redirect standard output and error

module load python3
export PYTHONPATH=$PYTHONPATH:/usr3/graduate/jaskielj/moments
readarray -t CMD_ARRAY < yft_bet_boot.modsel.runs
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub yft_bet_moments_array

#so everything get printed to yft_bet_boot.modsel

sed -i 's/\x00//g' yft_bet_boot.modsel
grep -a RESULT yft_bet_boot.modsel -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT  >yft_bet_boot.modsel.res1
#go through .res1 and make sure there aren't any lines that were cut off (delete those lines)
#Go through modSel_summary_yft_bet.R line by line; it doesn't really work if you run it as a job

NREPS=6 # number of random restarts per model per bootstrap rep
>winner.mods
for i in `seq 1 $NREPS`;do
awk '{print $2}' *.winmod | sed 's/$/.py/'  >>winner.mods;
done

sed -i -e 's/^/\/usr3\/graduate\/jaskielj\/AFS-analysis-with-moments\/work\//' winner.mods

CONTRAST=yft_bet_boot
ARGS="yft bet 154 91 0.0014 0.002"

>winner.modsel
for B in `seq 1 100`; do
INFILE=${CONTRAST}_${B}.sfs;
echo $INFILE;
NMODELS=`cat winner.mods | wc -l`
>winner.args
>${CONTRAST}.winboots
for i in `seq 1 $NMODELS`; do
echo "$INFILE $ARGS >>${CONTRAST}.winboots" >>winner.args;
done;
paste winner.mods winner.args -d " " >>winner.modsel;
done

#Running the winning model as an array
>yft_bet_winmodsel_array
nano yft_bet_winmodsel_array

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N yft_bet_winmodsel_array # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -pe omp 1
#$ -t 1-600
#$ -j y # Join standard output and error to a single file
#$ -o yft_bet_winmodsel_array.qlog # Name the file where to redirect standard output and error

module load python3
export PYTHONPATH=$PYTHONPATH:/usr3/graduate/jaskielj/moments
readarray -t CMD_ARRAY < winner.modsel
INDEX=$(($SGE_TASK_ID-1))
echo "Executing command: ${CMD_ARRAY[$INDEX]}"
eval "${CMD_ARRAY[$INDEX]}"

qsub yft_bet_winmodsel_array


## After Bootstrapping
#SKJ and YFT
CONTRAST=skj_yft_boot

#have to add -a to each grep because for some reason it reading as binary, I think because of weird ^@ characters?
#need to remove any lines that dont have all elements (just do in nano)
sed -i 's/\x00//g' ${CONTRAST}.winboots
grep -a RESULT ${CONTRAST}.winboots -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT >${CONTRAST}.winboots.res

##Run script bestBoot_summary_JJ.R to get parameter info
##See script demo_plots.R for plotting info
##FLAG: Add to main R file


#SKJ and BET
CONTRAST=skj_bet_boot

#have to add -a to each grep because for some reason it reading as binary, I think because of weird ^@ characters?
#need to remove any lines that dont have all elements (just do in nano)
sed -i 's/\x00//g' ${CONTRAST}.winboots
grep -a RESULT ${CONTRAST}.winboots -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT >${CONTRAST}.winboots.res

##Run script bestBoot_summary_JJ.R to get parameter info
##See script demo_plots.R for plotting info


#YFT and BET
CONTRAST=yft_bet_boot

#have to add -a to each grep because for some reason it reading as binary, I think because of weird ^@ characters?
#need to remove any lines that dont have all elements (just do in nano)
sed -i 's/\x00//g' ${CONTRAST}.winboots
grep -a RESULT ${CONTRAST}.winboots -A 4 | grep -aE "[0-9]|\]" | perl -pe 's/^100.+\.o\d+\S//' | perl -pe 's/\n//' | perl -pe 's/[\[\]]//g' | perl -pe 's/RESULT/\nRESULT/g' | grep -a RESULT >${CONTRAST}.winboots.res

##Run script bestBoot_summary_JJ.R to get parameter info
##See script demo_plots.R for plotting info

#-------------------------------------------Stairwayplot Analysis--------------------------------------------------
## Following the pipeline at https://github.com/xiaoming-liu/stairway-plot-v2
#download stairway_plot_v2.1.1.zip, add it to your home directory, and unzip
#unzipped files live in /usr3/graduate/jaskielj/Stairway_Plotting/stairway_plot_v2.1.1

## SFS work with individual species groups using LD thinned sites and folded (no anc or ref); well actually you need to unfold (include -anc), then it will fold the output

#Creating Folded SFS using the saf.idx files from before
realSFS skj_dem.saf.idx -fold 1 >skj_dem_folded.sfs
realSFS yft_dem.saf.idx -fold 1 >yft_dem_folded.sfs
realSFS bet_dem.saf.idx -fold 1 >bet_dem_folded.sfs

### Making blueprint files
#SKJ refgen w/o singletons
head -n1 skj_dem_folded.sfs | awk '{for(i=1;i<=203;i++) printf "%s ", $i; print ""}' > skj.folded.203.sfs #use this SFS in the blueprint

pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/skj_stairway

>skj.blueprint
nano skj.blueprint

#example blueprint file start
#input setting
popid: SKJ # id of the population (no white space)
nseq: 406 # number of sequences
L: 182061 # total number of observed nucleic sites, including polymorphic and monomorphic (get from running realSFS on SKJ.saf.idx without saving the output)
whether_folded: true # whethr the SFS is folded (true or false)
SFS: 167891.773875 6047.041004 1944.822591 1110.660164 714.765405 439.371985 482.937461 348.192693 136.712061 255.378353 207.329662 174.142639 203.642574 99.428730 52.544360 160.002358 109.082950 41.832165 36.427226 72.273843 68.910505 56.523909 60.735223 45.938574 36.865472 0.852424 38.850724 42.259062 37.696614 17.778467 5.826486 68.613986 29.381748 17.342846 26.829581 3.291455 0.163232 57.102365 10.592370 17.479765 6.467542 3.715302 36.344464 0.031415 0.032035 1.100105 41.482459 24.936243 0.119140 12.398015 4.274925 22.030500 1.518344 0.027910 14.393574 20.427504 0.040903 14.059095 20.341343 17.364596 2.233405 0.000091 0.000019 14.004799 0.365157 26.840915 0.053177 0.004087 0.003071 36.850901 0.796757 2.742814 13.067512 5.053696 0.060377 0.006225 0.117609 5.064116 25.551994 10.814951 1.202507 1.843513 19.190105 0.006408 0.000103 0.000010 0.003503 11.277350 1.029658 5.915666 15.330238 2.087506 0.254177 0.160605 0.756618 8.632986 17.290067 2.292569 0.221315 0.196600 0.989212 2.990516 28.427168 4.464256 0.006873 0.000198 0.002763 0.532653 9.341449 4.607843 2.561428 7.395518 21.158330 0.005052 0.033622 3.834682 16.472883 6.558659 4.663154 1.689079 0.229224 0.074669 0.165769 29.242282 10.345741 2.004231 0.037652 0.000567 0.003826 1.717913 14.903199 1.094961 0.417886 0.834094 2.711250 7.580618 13.453529 15.195739 10.266200 2.094143 0.051177 0.000188 0.000001 0.000000 0.000000 0.000002 0.007463 8.300504 15.362074 0.370005 0.303340 1.810059 11.423345 10.980332 0.278160 0.007787 0.003865 0.010261 0.027591 0.063156 0.383062 4.941962 4.579472 0.073202 0.002334 0.007922 0.245001 3.549430 13.873998 10.045600 2.248101 0.818724 0.794217 1.073297 0.814842 0.306668 0.248218 1.212683 6.909780 7.422017 2.199048 0.745583 1.180809 7.582818 6.126330 0.304199 0.156956 2.383939 12.956974 1.650308 0.021516 0.056758 10.353914 1.783798 0.523214 0.000382 0.032533 9.328013 9.240553 0.028893 0.000158 0.002658 2.146693 
smallest_size_of_SFS_bin_used_for_estimation: 2 # default is 1; to ignore singletons, uncomment this line and change this number to 2
largest_size_of_SFS_bin_used_for_estimation: 203 # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2
pct_training: 0.67 # percentage of sites for training
nrand: 10      50      100      200 # number of random break points for each try (separated by white space)
project_dir: /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/skj_stairway # project directory
stairway_plot_dir: /usr3/graduate/jaskielj/Stairway_Plotting/stairway_plot_v2.1.1/stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
random_seed: 8
#output setting
mu: 7.3e-9 # assumed mutation rate per site per generation
year_per_generation: 1 # assumed generation time (in years)
#plot setting
plot_title: SKJ # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,10000000 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
#example blueprint file end 

module load java
java -cp /usr3/graduate/jaskielj/Stairway_Plotting/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder skj.blueprint
#creates skj.blueprint.sh
nano skj.blueprint.sh #add header for job and submit

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N skj_blueprint # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 2
#$ -j y # Join standard output and error to a single file
#$ -o skj_blueprint.qlog # Name the file where to redirect standard output and error

qsub skj.blueprint.sh


#YFT
head -n1 yft_dem_folded.sfs | awk '{for(i=1;i<=96;i++) printf "%s ", $i; print ""}' > yft.folded.96.sfs

pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/yft_stairway

>yft.blueprint
nano yft.blueprint

#example blueprint file start
#input setting
popid: YFT # id of the population (no white space)
nseq: 193 # number of sequences
L: 158727 # total number of observed nucleic sites, including polymorphic and monomorphic (get from running realSFS on SKJ.saf.idx without saving the output)
whether_folded: true # whethr the SFS is folded (true or false)
SFS: 150595.151173 2953.700741 1484.511493 620.355047 626.320411 253.470211 261.737828 228.620039 99.407689 185.069519 41.133373 151.664651 43.960908 92.277699 45.306927 29.246956 64.421218 33.645176 50.190462 31.140652 16.060000 61.302696 9.738648 40.595215 84.080973 8.188220 23.093778 10.103289 9.902040 40.748805 39.082455 8.573391 24.094959 12.813090 44.568420 20.309394 15.343123 12.882461 21.011241 1.964373 10.952707 23.143967 28.682419 0.200471 0.029156 29.870008 0.000000 0.000003 6.756415 0.000000 0.000000 0.000000 4.803211 26.107175 7.409070 6.169610 0.000000 11.743526 0.625362 1.355271 16.810225 1.190843 0.000034 10.536838 7.823277 0.000055 0.000009 0.006812 14.124399 12.687138 0.000002 4.351056 0.000480 3.645965 3.576593 12.962081 0.445733 10.727187 5.271642 0.000822 0.000000 0.000000 0.000581 3.866145 13.093545 0.001691 21.275391 0.905194 7.008932 0.015556 1.036637 0.935427 1.254938 7.950402 11.857241 0.000014 
smallest_size_of_SFS_bin_used_for_estimation: 2 # default is 1; to ignore singletons, uncomment this line and change this number to 2
largest_size_of_SFS_bin_used_for_estimation: 96 # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2
pct_training: 0.67 # percentage of sites for training
nrand: 20     46      70      94 # number of random break points for each try (separated by white space)
project_dir: /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/yft_stairway # project directory
stairway_plot_dir: /usr3/graduate/jaskielj/Stairway_Plotting/stairway_plot_v2.1.1/stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
random_seed: 8
#output setting
mu: 7.3e-9 # assumed mutation rate per site per generation
year_per_generation: 2 # assumed generation time (in years)
#plot setting
plot_title: YFT # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,10000000 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
#example blueprint file end 

module load java
java -cp /usr3/graduate/jaskielj/Stairway_Plotting/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder yft.blueprint
#creates yft.blueprint.sh
nano yft.blueprint.sh #add header for job and submit

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N yft_blueprint # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 2
#$ -j y # Join standard output and error to a single file
#$ -o yft_blueprint.qlog # Name the file where to redirect standard output and error

qsub yft.blueprint.sh


#BET
pwd
/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/bet_stairway

#Needed to add this
head -n1 bet_dem_folded.sfs | awk '{for(i=1;i<=57;i++) printf "%s ", $i; print ""}' > bet.folded.57.sfs

>bet.blueprint
nano bet.blueprint

#example blueprint file start
#input setting
popid: BET # id of the population (no white space)
nseq: 115 # number of sequences
L: 167089 # total number of observed nucleic sites, including polymorphic and monomorphic (get from running realSFS on SKJ.saf.idx without saving the output)
whether_folded: true # whethr the SFS is folded (true or false)
SFS: 161743.561282 2610.444060 947.556464 454.594285 193.765567 226.012344 104.532161 62.665437 90.285445 64.253490 37.299410 33.759838 22.060784 21.350159 39.684073 18.654545 26.326794 18.842600 24.925505 16.484344 0.019564 23.617434 10.998683 8.675406 19.121469 4.797002 2.348036 4.863799 21.980257 11.843743 13.643037 0.002461 11.331491 1.948459 7.455264 19.724471 1.653598 18.685543 0.009950 24.126719 0.000000 10.153362 2.094187 0.213935 29.561056 0.000031 5.435324 7.809090 13.136315 8.211069 0.131222 22.513521 7.729778 6.480834 0.000001 0.000022 11.625277 
smallest_size_of_SFS_bin_used_for_estimation: 2 # default is 1; to ignore singletons, uncomment this line and change this number to 2
largest_size_of_SFS_bin_used_for_estimation: 57 # default is n-1; to ignore singletons, uncomment this line and change this number to nseq-2
pct_training: 0.67 # percentage of sites for training
nrand: 5      10      25      50 # number of random break points for each try (separated by white space)
project_dir: /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/bet_stairway # project directory
stairway_plot_dir: /usr3/graduate/jaskielj/Stairway_Plotting/stairway_plot_v2.1.1/stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
random_seed: 8
#output setting
mu: 7.3e-9 # assumed mutation rate per site per generation
year_per_generation: 3 # assumed generation time (in years)
#plot setting
plot_title: BET # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,10000000 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
#example blueprint file end 


module load java
java -cp /usr3/graduate/jaskielj/Stairway_Plotting/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder bet.blueprint
#creates yft.blueprint.sh
nano bet.blueprint.sh #add header for job and submit

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bet_blueprint # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 2
#$ -j y # Join standard output and error to a single file
#$ -o bet_blueprint.qlog # Name the file where to redirect standard output and error

qsub bet.blueprint.sh


##========================================= ANGSD Runs with Individual Species =======================================================
# Starting on these after removing low cov samples, removing TRs, and weird samples.

## - - - - - - - - - - - - - - - - - - - -  Just Skipjack - - - - - - - - - - - - - - - - - - -

# minind 162 (80% of 203) minQ 25 w/ -setMinDepthInd 5, LD thinned
>angsd_skj_mid5_noLD
nano angsd_skj_mid5_noLD

#!/bin/bash -l 
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_skj_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_skj_mid5_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -setMinDepthInd 5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 162 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBCF 1 -doPost 1 -doGlf 2"
angsd -b bams_skj -sites finalsites_0.2_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_skj_mid5_noLD

qsub angsd_skj_mid5_noLD

# how many SNPs?
NSITES=`zcat angsd_skj_mid5_noLD.mafs.gz | wc -l`
echo $NSITES
#334


## - - - - - - - - - - - - - - - - - - - -  Just Yellowfin - - - - - - - - - - - - - - - - - - -

# minind 77 (80% of 96) minQ 25 w/ -setMinDepthInd 5, LD thinned
>angsd_yft_mid5_noLD
nano angsd_yft_mid5_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_yft_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_yft_mid5_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 77 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_yft -sites finalsites_0.2_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_yft_mid5_noLD

qsub angsd_yft_mid5_noLD

NSITES=`zcat angsd_yft_mid5_noLD.mafs.gz | wc -l` 
echo $NSITES
#383


# - - - - - - - - - - - - - - - - - - - - Just Bigeye - - - - - - - - - - - - - - - - - - - - - - - 

# minind 45 (80% of 57) minQ 25 w/ -setMinDepthInd 5, LD thinned
>angsd_bet_mid5_noLD
nano angsd_bet_mid5_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_bet_mid5_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_bet_mid5_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 45 -setMinDepthInd 5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2" 
angsd -b bams_bet -sites finalsites_0.2_unlinked -GL 1 $FILTERS $TODO -P 1 -out angsd_bet_mid5_noLD

qsub angsd_bet_mid5_noLD

NSITES=`zcat angsd_bet_mid5_noLD.mafs.gz | wc -l` 
echo $NSITES
#302


# ======================================= Basic Summary Stats ======================================

# - - - - - - - - - - - - - - - - - - - - - - - Heterozygosity - - - - - - - - - - - - - - - - - - - - 
# Skipjack (203 individuals, 80% is 162), LD thinned sites
>angsd_het_skj_noLD
nano angsd_het_skj_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_skj_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_skj_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 162"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_skj -sites finalsites_0.2_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_skj_noLD

qsub angsd_het_skj_noLD

# scp .beagle.gz and use R script to calculate heterozygosity and plot


# Yellowfin (96 individuals, 80% is 77), LD thinned sites
>angsd_het_yft_noLD
nano angsd_het_yft_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_yft_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_yft_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 77"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_yft -sites finalsites_0.2_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_yft_noLD

qsub angsd_het_yft_noLD


# Bigeye (57 individuals, 80% is 46) LD thinned
>angsd_het_bet_noLD
nano angsd_het_bet_noLD

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_het_bet_noLD # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o angsd_het_bet_noLD.qlog # Name the file where to redirect standard output and error

module load angsd
FILTERS="-maxHetFreq 0.5 -uniqueOnly 1 -remove_bads 1 -skipTriallelic 1 -minMapQ 20 -minQ 25 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -minInd 46"
TODO="-makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 32 -doPost 1 -doGlf 2 -doCounts 1 -doMajorMinor 1 -dosnpstat 1 -doMaf 1"
angsd -bam bams_bet -sites finalsites_0.2_unlinked -GL 1 -P 1 $TODO $FILTERS -out angsd_het_bet_noLD

qsub angsd_het_bet_noLD


#------------------------------- Pairwise Fst between years of sampling within species --------------------------------------
#first need to make individual BAM lists for each year of sampling for each species

#Skipjack (ref)
#Doing SAF/SFS work

#2015 (17 individuals, so 80% is 14)
>sfs_skj_2015
nano sfs_skj_2015

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_skj_2015 -GL 1 -P 1 -minInd 14 $TODO -out skj_2015_allsites

qsub sfs_skj_2015


#2016 (69 individuals so 80% is 55)
>sfs_skj_2016
nano sfs_skj_2016

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_skj_2016 -GL 1 -P 1 -minInd 55 $TODO -out skj_2016_allsites

qsub sfs_skj_2016


#2017 (27 individuals so 80% is 22)
>sfs_skj_2017
nano sfs_skj_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017.qlog # Name the file where to redirect standard output and error

module load angsd

export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_skj_2017 -GL 1 -P 1 -minInd 22 $TODO -out skj_2017_allsites

qsub sfs_skj_2017


#2018 (73 individuals so 80% is 58)
>sfs_skj_2018
nano sfs_skj_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2018.qlog # Name the file where to redirect standard output and error

module load angsd

export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_skj_2018 -GL 1 -P 1 -minInd 58 $TODO -out skj_2018_allsites

qsub sfs_skj_2018


#2019 (16 individuals so 80% is 13)
>sfs_skj_2019
nano sfs_skj_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2019.qlog # Name the file where to redirect standard output and error

module load angsd

export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_skj_2019 -GL 1 -P 1 -minInd 13 $TODO -out skj_2019_allsites

qsub sfs_skj_2019


#2022 (only 1 individual so no minind)
>sfs_skj_2022
nano sfs_skj_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2022.qlog # Name the file where to redirect standard output and error

module load angsd

export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_skj_2022 -GL 1 -P 1 -minInd 0 $TODO -out skj_2022_allsites

qsub sfs_skj_2022


# Generating per-population SFS
realSFS skj_2015_allsites.saf.idx >skj_2015_allsites.sfs
realSFS skj_2016_allsites.saf.idx >skj_2016_allsites.sfs
realSFS skj_2017_allsites.saf.idx >skj_2017_allsites.sfs
realSFS skj_2018_allsites.saf.idx >skj_2018_allsites.sfs
realSFS skj_2019_allsites.saf.idx >skj_2019_allsites.sfs
realSFS skj_2022_allsites.saf.idx >skj_2022_allsites.sfs

# Writing down 2d-SFS priors - 2015 & 2016
>sfs_skj_2015_2016
nano sfs_skj_2015_2016

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2016 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2016.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2015_allsites.saf.idx skj_2016_allsites.saf.idx -P 1 > skj_2015_2016_allsites.sfs ; realSFS fst index skj_2015_allsites.saf.idx skj_2016_allsites.saf.idx -sfs skj_2015_2016_allsites.sfs -fstout skj_2015_2016_allsites

qsub sfs_skj_2015_2016

# Global Fst between populations 
realSFS fst stats skj_2015_2016_allsites.fst.idx

# Output:
FST.Unweight[nObs:179810]:0.002991 Fst.Weight:0.011931

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# Writing down 2d-SFS priors - 2015 & 2017
>sfs_skj_2015_2017
nano sfs_skj_2015_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2017.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2015_allsites.saf.idx skj_2017_allsites.saf.idx -P 1 > skj_2015_2017_allsites.sfs ; realSFS fst index skj_2015_allsites.saf.idx skj_2017_allsites.saf.idx -sfs skj_2015_2017_allsites.sfs -fstout skj_2015_2017_allsites

qsub sfs_skj_2015_2017

# Global Fst between populations 
realSFS fst stats skj_2015_2017_allsites.fst.idx

# Output:
FST.Unweight[nObs:177772]:0.008951 Fst.Weight:0.017096


# Writing down 2d-SFS priors - 2015 & 2018
>sfs_skj_2015_2018
nano sfs_skj_2015_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2015_allsites.saf.idx skj_2018_allsites.saf.idx -P 1 > skj_2015_2018_allsites.sfs ; realSFS fst index skj_2015_allsites.saf.idx skj_2018_allsites.saf.idx -sfs skj_2015_2018_allsites.sfs -fstout skj_2015_2018_allsites

qsub sfs_skj_2015_2018

# Global Fst between populations 
realSFS fst stats skj_2015_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:182917]:0.003578 Fst.Weight:0.009042


# Writing down 2d-SFS priors - 2015 & 2019
>sfs_skj_2015_2019
nano sfs_skj_2015_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2015_allsites.saf.idx skj_2019_allsites.saf.idx -P 1 > skj_2015_2019_allsites.sfs ; realSFS fst index skj_2015_allsites.saf.idx skj_2019_allsites.saf.idx -sfs skj_2015_2019_allsites.sfs -fstout skj_2015_2019_allsites

qsub sfs_skj_2015_2019

# Global Fst between populations 
realSFS fst stats skj_2015_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:91671]:0.018474 Fst.Weight:0.020767


# Writing down 2d-SFS priors - 2015 & 2022
>sfs_skj_2015_2022
nano sfs_skj_2015_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2015_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2015_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2015_allsites.saf.idx skj_2022_allsites.saf.idx -P 1 > skj_2015_2022_allsites.sfs ; realSFS fst index skj_2015_allsites.saf.idx skj_2022_allsites.saf.idx -sfs skj_2015_2022_allsites.sfs -fstout skj_2015_2022_allsites

qsub sfs_skj_2015_2022

# Global Fst between populations 
realSFS fst stats skj_2015_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:186051]:0.061446 Fst.Weight:0.207243


# Writing down 2d-SFS priors - 2016 & 2017
>sfs_skj_2016_2017
nano sfs_skj_2016_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2017.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2016_allsites.saf.idx skj_2017_allsites.saf.idx -P 1 > skj_2016_2017_allsites.sfs ; realSFS fst index skj_2016_allsites.saf.idx skj_2017_allsites.saf.idx -sfs skj_2016_2017_allsites.sfs -fstout skj_2016_2017_allsites

qsub sfs_skj_2016_2017

# Global Fst between populations 
realSFS fst stats skj_2016_2017_allsites.fst.idx

# Output:
FST.Unweight[nObs:175963]:0.006348 Fst.Weight:0.006625


# Writing down 2d-SFS priors - 2016 & 2018
>sfs_skj_2016_2018
nano sfs_skj_2016_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2016_allsites.saf.idx skj_2018_allsites.saf.idx -P 1 > skj_2016_2018_allsites.sfs ; realSFS fst index skj_2016_allsites.saf.idx skj_2018_allsites.saf.idx -sfs skj_2016_2018_allsites.sfs -fstout skj_2016_2018_allsites

qsub sfs_skj_2016_2018

# Global Fst between populations 
realSFS fst stats skj_2016_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:180439]:0.003599 Fst.Weight:0.004088


# Writing down 2d-SFS priors - 2016 & 2019
>sfs_skj_2016_2019
nano sfs_skj_2016_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2016_allsites.saf.idx skj_2019_allsites.saf.idx -P 1 > skj_2016_2019_allsites.sfs ; realSFS fst index skj_2016_allsites.saf.idx skj_2019_allsites.saf.idx -sfs skj_2016_2019_allsites.sfs -fstout skj_2016_2019_allsites

qsub sfs_skj_2016_2019

# Global Fst between populations 
realSFS fst stats skj_2016_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:91687]:0.023775 Fst.Weight:0.013059


# Writing down 2d-SFS priors - 2016 & 2022
>sfs_skj_2016_2022
nano sfs_skj_2016_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2016_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2016_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2016_allsites.saf.idx skj_2022_allsites.saf.idx -P 1 > skj_2016_2022_allsites.sfs ; realSFS fst index skj_2016_allsites.saf.idx skj_2022_allsites.saf.idx -sfs skj_2016_2022_allsites.sfs -fstout skj_2016_2022_allsites

qsub sfs_skj_2016_2022

# Global Fst between populations 
realSFS fst stats skj_2016_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:180123]:-0.045005 Fst.Weight:0.239060


# Writing down 2d-SFS priors - 2017 & 2018
>sfs_skj_2017_2018
nano sfs_skj_2017_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2017_allsites.saf.idx skj_2018_allsites.saf.idx -P 1 > skj_2017_2018_allsites.sfs ; realSFS fst index skj_2017_allsites.saf.idx skj_2018_allsites.saf.idx -sfs skj_2017_2018_allsites.sfs -fstout skj_2017_2018_allsites

qsub sfs_skj_2017_2018

# Global Fst between populations 
realSFS fst stats skj_2017_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:178535]:0.007089 Fst.Weight:0.008180


# Writing down 2d-SFS priors - 2017 & 2019
>sfs_skj_2017_2019
nano sfs_skj_2017_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2017_allsites.saf.idx skj_2019_allsites.saf.idx -P 1 > skj_2017_2019_allsites.sfs ; realSFS fst index skj_2017_allsites.saf.idx skj_2019_allsites.saf.idx -sfs skj_2017_2019_allsites.sfs -fstout skj_2017_2019_allsites

qsub sfs_skj_2017_2019

# Global Fst between populations 
realSFS fst stats skj_2017_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:91413]:0.020455 Fst.Weight:0.015080


# Writing down 2d-SFS priors - 2017 & 2022
>sfs_skj_2017_2022
nano sfs_skj_2017_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2017_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2017_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2017_allsites.saf.idx skj_2022_allsites.saf.idx -P 1 > skj_2017_2022_allsites.sfs ; realSFS fst index skj_2017_allsites.saf.idx skj_2022_allsites.saf.idx -sfs skj_2017_2022_allsites.sfs -fstout skj_2017_2022_allsites

qsub sfs_skj_2017_2022

# Global Fst between populations 
realSFS fst stats skj_2017_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:178508]:-0.029137 Fst.Weight:0.257691


# Writing down 2d-SFS priors - 2018 & 2019
>sfs_skj_2018_2019
nano sfs_skj_2018_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2018_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2018_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2018_allsites.saf.idx skj_2019_allsites.saf.idx -P 1 > skj_2018_2019_allsites.sfs ; realSFS fst index skj_2018_allsites.saf.idx skj_2019_allsites.saf.idx -sfs skj_2018_2019_allsites.sfs -fstout skj_2018_2019_allsites

qsub sfs_skj_2018_2019

# Global Fst between populations 
realSFS fst stats skj_2018_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:91862]:0.025652 Fst.Weight:0.013550


# Writing down 2d-SFS priors - 2018 & 2022
>sfs_skj_2018_2022
nano sfs_skj_2018_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2018_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2018_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2018_allsites.saf.idx skj_2022_allsites.saf.idx -P 1 > skj_2018_2022_allsites.sfs ; realSFS fst index skj_2018_allsites.saf.idx skj_2022_allsites.saf.idx -sfs skj_2018_2022_allsites.sfs -fstout skj_2018_2022_allsites

qsub sfs_skj_2018_2022

# Global Fst between populations 
realSFS fst stats skj_2018_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:183959]:-0.043227 Fst.Weight:0.218430


# Writing down 2d-SFS priors - 2019 & 2022
>sfs_skj_2019_2022
nano sfs_skj_2019_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_skj_2019_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_skj_2019_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS skj_2019_allsites.saf.idx skj_2022_allsites.saf.idx -P 1 > skj_2019_2022_allsites.sfs ; realSFS fst index skj_2019_allsites.saf.idx skj_2022_allsites.saf.idx -sfs skj_2019_2022_allsites.sfs -fstout skj_2019_2022_allsites

qsub sfs_skj_2019_2022

# Global Fst between populations 
realSFS fst stats skj_2019_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:91677]:-0.079070 Fst.Weight:0.244375


#------------------------ Yellowfin reference-based analysis of genetic divergence (FST) across years
cd /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

#Doing SAF/SFS work

#2015 (8 individuals, so 80% is 6)
>sfs_yft_2015
nano sfs_yft_2015

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft_2015 -GL 1 -P 1 -minInd 6 $TODO -out yft_2015_allsites

qsub sfs_yft_2015


#2016 (4 individuals, so 80% is 3)
>sfs_yft_2016
nano sfs_yft_2016

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft_2016 -GL 1 -P 1 -minInd 3 $TODO -out yft_2016_allsites

qsub sfs_yft_2016


#2017 (5 individuals, so 80% is 4)
>sfs_yft_2017
nano sfs_yft_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft_2017 -GL 1 -P 1 -minInd 4 $TODO -out yft_2017_allsites

qsub sfs_yft_2017


#2018 (43 individuals, so 80% is 34)
>sfs_yft_2018
nano sfs_yft_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft_2018 -GL 1 -P 1 -minInd 34 $TODO -out yft_2018_allsites

qsub sfs_yft_2018


#2019 (1 individual, so no minind)
>sfs_yft_2019
nano sfs_yft_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2019.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft_2019 -GL 1 -P 1 -minInd 0 $TODO -out yft_2019_allsites

qsub sfs_yft_2019


#2020 (29 individuals, so 80% is 23)
>sfs_yft_2020
nano sfs_yft_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2020.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft_2020 -GL 1 -P 1 -minInd 23 $TODO -out yft_2020_allsites

qsub sfs_yft_2020


#2022 (6 individuals, so 80% is 5)
>sfs_yft_2022
nano sfs_yft_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2022.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_yft_2022 -GL 1 -P 1 -minInd 5 $TODO -out yft_2022_allsites

qsub sfs_yft_2022


# Generating per-population SFS
realSFS yft_2015_allsites.saf.idx >yft_2015_allsites.sfs
realSFS yft_2016_allsites.saf.idx >yft_2016_allsites.sfs
realSFS yft_2017_allsites.saf.idx >yft_2017_allsites.sfs
realSFS yft_2018_allsites.saf.idx >yft_2018_allsites.sfs
realSFS yft_2019_allsites.saf.idx >yft_2019_allsites.sfs
realSFS yft_2020_allsites.saf.idx >yft_2020_allsites.sfs
realSFS yft_2022_allsites.saf.idx >yft_2022_allsites.sfs

# Writing down 2d-SFS priors - 2015 & 2016
>sfs_yft_2015_2016
nano sfs_yft_2015_2016

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2016 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2016.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2015_allsites.saf.idx yft_2016_allsites.saf.idx -P 1 > yft_2015_2016_allsites.sfs ; realSFS fst index yft_2015_allsites.saf.idx yft_2016_allsites.saf.idx -sfs yft_2015_2016_allsites.sfs -fstout yft_2015_2016_allsites

qsub sfs_yft_2015_2016

# Global Fst between populations 
realSFS fst stats yft_2015_2016_allsites.fst.idx

# Output:
FST.Unweight[nObs:137647]:0.027010 Fst.Weight:0.057093

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# Writing down 2d-SFS priors - 2015 & 2017
>sfs_yft_2015_2017
nano sfs_yft_2015_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2017.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2015_allsites.saf.idx yft_2017_allsites.saf.idx -P 1 > yft_2015_2017_allsites.sfs ; realSFS fst index yft_2015_allsites.saf.idx yft_2017_allsites.saf.idx -sfs yft_2015_2017_allsites.sfs -fstout yft_2015_2017_allsites

qsub sfs_yft_2015_2017

# Global Fst between populations 
realSFS fst stats yft_2015_2017_allsites.fst.idx

# Output:
FST.Unweight[nObs:99014]:0.056400 Fst.Weight:0.061763


# Writing down 2d-SFS priors - 2015 & 2018
>sfs_yft_2015_2018
nano sfs_yft_2015_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2015_allsites.saf.idx yft_2018_allsites.saf.idx -P 1 > yft_2015_2018_allsites.sfs ; realSFS fst index yft_2015_allsites.saf.idx yft_2018_allsites.saf.idx -sfs yft_2015_2018_allsites.sfs -fstout yft_2015_2018_allsites

qsub sfs_yft_2015_2018

# Global Fst between populations 
realSFS fst stats yft_2015_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:136123]:0.040299 Fst.Weight:0.028922


# Writing down 2d-SFS priors - 2015 & 2019
>sfs_yft_2015_2019
nano sfs_yft_2015_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2015_allsites.saf.idx yft_2019_allsites.saf.idx -P 1 > yft_2015_2019_allsites.sfs ; realSFS fst index yft_2015_allsites.saf.idx yft_2019_allsites.saf.idx -sfs yft_2015_2019_allsites.sfs -fstout yft_2015_2019_allsites

qsub sfs_yft_2015_2019

# Global Fst between populations 
realSFS fst stats yft_2015_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:114212]:0.080622 Fst.Weight:0.172294


# Writing down 2d-SFS priors - 2015 & 2020
>sfs_yft_2015_2020
nano sfs_yft_2015_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2015_allsites.saf.idx yft_2020_allsites.saf.idx -P 1 > yft_2015_2020_allsites.sfs ; realSFS fst index yft_2015_allsites.saf.idx yft_2020_allsites.saf.idx -sfs yft_2015_2020_allsites.sfs -fstout yft_2015_2020_allsites

qsub sfs_yft_2015_2020

# Global Fst between populations 
realSFS fst stats yft_2015_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:127740]:0.038660 Fst.Weight:0.028747


# Writing down 2d-SFS priors - 2015 & 2022
>sfs_yft_2015_2022
nano sfs_yft_2015_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2015_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2015_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2015_allsites.saf.idx yft_2022_allsites.saf.idx -P 1 > yft_2015_2022_allsites.sfs ; realSFS fst index yft_2015_allsites.saf.idx yft_2022_allsites.saf.idx -sfs yft_2015_2022_allsites.sfs -fstout yft_2015_2022_allsites

qsub sfs_yft_2015_2022

# Global Fst between populations 
realSFS fst stats yft_2015_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:140858]:0.025742 Fst.Weight:0.050722


# Writing down 2d-SFS priors - 2016 & 2017
>sfs_yft_2016_2017
nano sfs_yft_2016_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2017.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2016_allsites.saf.idx yft_2017_allsites.saf.idx -P 1 > yft_2016_2017_allsites.sfs ; realSFS fst index yft_2016_allsites.saf.idx yft_2017_allsites.saf.idx -sfs yft_2016_2017_allsites.sfs -fstout yft_2016_2017_allsites

qsub sfs_yft_2016_2017

# Global Fst between populations 
realSFS fst stats yft_2016_2017_allsites.fst.idx

# Output:
FST.Unweight[nObs:115820]:0.065605 Fst.Weight:0.067813


# Writing down 2d-SFS priors - 2016 & 2018
>sfs_yft_2016_2018
nano sfs_yft_2016_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2016_allsites.saf.idx yft_2018_allsites.saf.idx -P 1 > yft_2016_2018_allsites.sfs ; realSFS fst index yft_2016_allsites.saf.idx yft_2018_allsites.saf.idx -sfs yft_2016_2018_allsites.sfs -fstout yft_2016_2018_allsites

qsub sfs_yft_2016_2018

# Global Fst between populations 
realSFS fst stats yft_2016_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:171060]:0.025108 Fst.Weight:0.034085


# Writing down 2d-SFS priors - 2016 & 2019
>sfs_yft_2016_2019
nano sfs_yft_2016_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2016_allsites.saf.idx yft_2019_allsites.saf.idx -P 1 > yft_2016_2019_allsites.sfs ; realSFS fst index yft_2016_allsites.saf.idx yft_2019_allsites.saf.idx -sfs yft_2016_2019_allsites.sfs -fstout yft_2016_2019_allsites

qsub sfs_yft_2016_2019

# Global Fst between populations 
realSFS fst stats yft_2016_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:142388]:0.281500 Fst.Weight:0.211484


# Writing down 2d-SFS priors - 2016 & 2020
>sfs_yft_2016_2020
nano sfs_yft_2016_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2016_allsites.saf.idx yft_2020_allsites.saf.idx -P 1 > yft_2016_2020_allsites.sfs ; realSFS fst index yft_2016_allsites.saf.idx yft_2020_allsites.saf.idx -sfs yft_2016_2020_allsites.sfs -fstout yft_2016_2020_allsites

qsub sfs_yft_2016_2020

# Global Fst between populations 
realSFS fst stats yft_2016_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:158880]:0.033125 Fst.Weight:0.039679


# Writing down 2d-SFS priors - 2016 & 2022
>sfs_yft_2016_2022
nano sfs_yft_2016_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2016_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2016_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2016_allsites.saf.idx yft_2022_allsites.saf.idx -P 1 > yft_2016_2022_allsites.sfs ; realSFS fst index yft_2016_allsites.saf.idx yft_2022_allsites.saf.idx -sfs yft_2016_2022_allsites.sfs -fstout yft_2016_2022_allsites

qsub sfs_yft_2016_2022

# Global Fst between populations 
realSFS fst stats yft_2016_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:184853]:0.091526 Fst.Weight:0.055142


# Writing down 2d-SFS priors - 2017 & 2018
>sfs_yft_2017_2018
nano sfs_yft_2017_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2017_allsites.saf.idx yft_2018_allsites.saf.idx -P 1 > yft_2017_2018_allsites.sfs ; realSFS fst index yft_2017_allsites.saf.idx yft_2018_allsites.saf.idx -sfs yft_2017_2018_allsites.sfs -fstout yft_2017_2018_allsites

qsub sfs_yft_2017_2018

# Global Fst between populations 
realSFS fst stats yft_2017_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:114628]:0.098369 Fst.Weight:0.045528


# Writing down 2d-SFS priors - 2017 & 2019
>sfs_yft_2017_2019
nano sfs_yft_2017_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2017_allsites.saf.idx yft_2019_allsites.saf.idx -P 1 > yft_2017_2019_allsites.sfs ; realSFS fst index yft_2017_allsites.saf.idx yft_2019_allsites.saf.idx -sfs yft_2017_2019_allsites.sfs -fstout yft_2017_2019_allsites

qsub sfs_yft_2017_2019

# Global Fst between populations 
realSFS fst stats yft_2017_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:99002]:0.131173 Fst.Weight:0.160696


# Writing down 2d-SFS priors - 2017 & 2020
>sfs_yft_2017_2020
nano sfs_yft_2017_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2017_allsites.saf.idx yft_2020_allsites.saf.idx -P 1 > yft_2017_2020_allsites.sfs ; realSFS fst index yft_2017_allsites.saf.idx yft_2020_allsites.saf.idx -sfs yft_2017_2020_allsites.sfs -fstout yft_2017_2020_allsites

qsub sfs_yft_2017_2020

# Global Fst between populations 
realSFS fst stats yft_2017_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:108388]:0.122381 Fst.Weight:0.057365


# Writing down 2d-SFS priors - 2017 & 2022
>sfs_yft_2017_2022
nano sfs_yft_2017_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2017_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2017_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2017_allsites.saf.idx yft_2022_allsites.saf.idx -P 1 > yft_2017_2022_allsites.sfs ; realSFS fst index yft_2017_allsites.saf.idx yft_2022_allsites.saf.idx -sfs yft_2017_2022_allsites.sfs -fstout yft_2017_2022_allsites

qsub sfs_yft_2017_2022

# Global Fst between populations 
realSFS fst stats yft_2017_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:118507]:0.101022 Fst.Weight:0.054234


# Writing down 2d-SFS priors - 2018 & 2019
>sfs_yft_2018_2019
nano sfs_yft_2018_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2018_allsites.saf.idx yft_2019_allsites.saf.idx -P 1 > yft_2018_2019_allsites.sfs ; realSFS fst index yft_2018_allsites.saf.idx yft_2019_allsites.saf.idx -sfs yft_2018_2019_allsites.sfs -fstout yft_2018_2019_allsites

qsub sfs_yft_2018_2019

# Global Fst between populations 
realSFS fst stats yft_2018_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:138835]:0.177493 Fst.Weight:0.145349


# Writing down 2d-SFS priors - 2018 & 2020
>sfs_yft_2018_2020
nano sfs_yft_2018_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2018_allsites.saf.idx yft_2020_allsites.saf.idx -P 1 > yft_2018_2020_allsites.sfs ; realSFS fst index yft_2018_allsites.saf.idx yft_2020_allsites.saf.idx -sfs yft_2018_2020_allsites.sfs -fstout yft_2018_2020_allsites

qsub sfs_yft_2018_2020

# Global Fst between populations 
realSFS fst stats yft_2018_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:156893]:0.009547 Fst.Weight:0.014405


# Writing down 2d-SFS priors - 2018 & 2022
>sfs_yft_2018_2022
nano sfs_yft_2018_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2018_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2018_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2018_allsites.saf.idx yft_2022_allsites.saf.idx -P 1 > yft_2018_2022_allsites.sfs ; realSFS fst index yft_2018_allsites.saf.idx yft_2022_allsites.saf.idx -sfs yft_2018_2022_allsites.sfs -fstout yft_2018_2022_allsites

qsub sfs_yft_2018_2022

# Global Fst between populations 
realSFS fst stats yft_2018_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:176964]:-0.005932 Fst.Weight:0.026472


# Writing down 2d-SFS priors - 2019 & 2020
>sfs_yft_2019_2020
nano sfs_yft_2019_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2019_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2019_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2019_allsites.saf.idx yft_2020_allsites.saf.idx -P 1 > yft_2019_2020_allsites.sfs ; realSFS fst index yft_2019_allsites.saf.idx yft_2020_allsites.saf.idx -sfs yft_2019_2020_allsites.sfs -fstout yft_2019_2020_allsites

qsub sfs_yft_2019_2020

# Global Fst between populations 
realSFS fst stats yft_2019_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:130176]:0.150524 Fst.Weight:0.148192


# Writing down 2d-SFS priors - 2019 & 2022
>sfs_yft_2019_2022
nano sfs_yft_2019_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2019_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2019_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2019_allsites.saf.idx yft_2022_allsites.saf.idx -P 1 > yft_2019_2022_allsites.sfs ; realSFS fst index yft_2019_allsites.saf.idx yft_2022_allsites.saf.idx -sfs yft_2019_2022_allsites.sfs -fstout yft_2019_2022_allsites

qsub sfs_yft_2019_2022

# Global Fst between populations 
realSFS fst stats yft_2019_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:146870]:0.612832 Fst.Weight:0.173028


# Writing down 2d-SFS priors - 2020 & 2022
>sfs_yft_2020_2022
nano sfs_yft_2020_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_yft_2020_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_yft_2020_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS yft_2020_allsites.saf.idx yft_2022_allsites.saf.idx -P 1 > yft_2020_2022_allsites.sfs ; realSFS fst index yft_2020_allsites.saf.idx yft_2022_allsites.saf.idx -sfs yft_2020_2022_allsites.sfs -fstout yft_2020_2022_allsites

qsub sfs_yft_2020_2022

# Global Fst between populations 
realSFS fst stats yft_2020_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:163753]:-0.003248 Fst.Weight:0.037041


#------------------------ Bigeye reference-based analysis of genetic divergence (FST) across years
cd /projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis

#Doing SAF/SFS work

#2016 (11 individuals, so 80% is 9)
>sfs_bet_2016
nano sfs_bet_2016

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_bet_2016 -GL 1 -P 1 -minInd 9 $TODO -out bet_2016_allsites

qsub sfs_bet_2016


#2017 (2 individuals, so set minind to 1)
>sfs_bet_2017
nano sfs_bet_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_bet_2017 -GL 1 -P 1 -minInd 1 $TODO -out bet_2017_allsites

qsub sfs_bet_2017


#2018 (11 individuals, so set minind to 9)
>sfs_bet_2018
nano sfs_bet_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_bet_2018 -GL 1 -P 1 -minInd 9 $TODO -out bet_2018_allsites

qsub sfs_bet_2018


#2019 (1 individual, so set minind to 0)
>sfs_bet_2019
nano sfs_bet_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2019.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_bet_2019 -GL 1 -P 1 -minInd 0 $TODO -out bet_2019_allsites

qsub sfs_bet_2019


#2020 (24 individuals, so set minind to 19)
>sfs_bet_2020
nano sfs_bet_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2020.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_bet_2020 -GL 1 -P 1 -minInd 19 $TODO -out bet_2020_allsites

qsub sfs_bet_2020


#2022 (8 individuals, so set minind to 6)
>sfs_bet_2022
nano sfs_bet_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2022.qlog # Name the file where to redirect standard output and error

module load angsd
export GENOME_REF=GCF_914725855.1_fThuAlb1.1_genomic.fasta
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF"
angsd -sites finalsites_0.2_unlinked -b bams_bet_2022 -GL 1 -P 1 -minInd 6 $TODO -out bet_2022_allsites

qsub sfs_bet_2022

# Generating per-population SFS
realSFS bet_2016_allsites.saf.idx >bet_2016_allsites.sfs
realSFS bet_2017_allsites.saf.idx >bet_2017_allsites.sfs
realSFS bet_2018_allsites.saf.idx >bet_2018_allsites.sfs
realSFS bet_2019_allsites.saf.idx >bet_2019_allsites.sfs
realSFS bet_2020_allsites.saf.idx >bet_2020_allsites.sfs
realSFS bet_2022_allsites.saf.idx >bet_2022_allsites.sfs

# Writing down 2d-SFS priors - 2016 & 2017
>sfs_bet_2016_2017
nano sfs_bet_2016_2017

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2017 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2017.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2016_allsites.saf.idx bet_2017_allsites.saf.idx -P 1 > bet_2016_2017_allsites.sfs ; realSFS fst index bet_2016_allsites.saf.idx bet_2017_allsites.saf.idx -sfs bet_2016_2017_allsites.sfs -fstout bet_2016_2017_allsites

qsub sfs_bet_2016_2017

# Global Fst between populations 
realSFS fst stats bet_2016_2017_allsites.fst.idx

# Output:
FST.Unweight[nObs:125382]:0.004797 Fst.Weight:0.074624

# The difference between unweighted and weighted values is averaging versus ratio of sums method. 
# The weighted value is derived from the ratio of the separate per-locus sums of numerator and denominator values, 
# while the unweighted value is the average of per-locus values. [If that is unclear: weighted is sum(a)/sum(a+b), while unweighted is average(a/(a+b))].


# Writing down 2d-SFS priors - 2016 & 2018
>sfs_bet_2016_2018
nano sfs_bet_2016_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2016_allsites.saf.idx bet_2018_allsites.saf.idx -P 1 > bet_2016_2018_allsites.sfs ; realSFS fst index bet_2016_allsites.saf.idx bet_2018_allsites.saf.idx -sfs bet_2016_2018_allsites.sfs -fstout bet_2016_2018_allsites

qsub sfs_bet_2016_2018

# Global Fst between populations 
realSFS fst stats bet_2016_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:125229]:0.024581 Fst.Weight:0.029608


# Writing down 2d-SFS priors - 2016 & 2019
>sfs_bet_2016_2019
nano sfs_bet_2016_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2016_allsites.saf.idx bet_2019_allsites.saf.idx -P 1 > bet_2016_2019_allsites.sfs ; realSFS fst index bet_2016_allsites.saf.idx bet_2019_allsites.saf.idx -sfs bet_2016_2019_allsites.sfs -fstout bet_2016_2019_allsites

qsub sfs_bet_2016_2019

# Global Fst between populations 
realSFS fst stats bet_2016_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:76848]:0.253198 Fst.Weight:0.173878


# Writing down 2d-SFS priors - 2016 & 2020
>sfs_bet_2016_2020
nano sfs_bet_2016_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2016_allsites.saf.idx bet_2020_allsites.saf.idx -P 1 > bet_2016_2020_allsites.sfs ; realSFS fst index bet_2016_allsites.saf.idx bet_2020_allsites.saf.idx -sfs bet_2016_2020_allsites.sfs -fstout bet_2016_2020_allsites

qsub sfs_bet_2016_2020

# Global Fst between populations 
realSFS fst stats bet_2016_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:115059]:0.021960 Fst.Weight:0.023491


# Writing down 2d-SFS priors - 2016 & 2022
>sfs_bet_2016_2022
nano sfs_bet_2016_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2016_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2016_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2016_allsites.saf.idx bet_2022_allsites.saf.idx -P 1 > bet_2016_2022_allsites.sfs ; realSFS fst index bet_2016_allsites.saf.idx bet_2022_allsites.saf.idx -sfs bet_2016_2022_allsites.sfs -fstout bet_2016_2022_allsites

qsub sfs_bet_2016_2022

# Global Fst between populations 
realSFS fst stats bet_2016_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:125647]:0.016662 Fst.Weight:0.034225


# Writing down 2d-SFS priors - 2017 & 2018
>sfs_bet_2017_2018
nano sfs_bet_2017_2018

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2018 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2018.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2017_allsites.saf.idx bet_2018_allsites.saf.idx -P 1 > bet_2017_2018_allsites.sfs ; realSFS fst index bet_2017_allsites.saf.idx bet_2018_allsites.saf.idx -sfs bet_2017_2018_allsites.sfs -fstout bet_2017_2018_allsites

qsub sfs_bet_2017_2018

# Global Fst between populations 
realSFS fst stats bet_2017_2018_allsites.fst.idx

# Output:
FST.Unweight[nObs:187526]:0.115034 Fst.Weight:0.079525


# Writing down 2d-SFS priors - 2017 & 2019
>sfs_bet_2017_2019
nano sfs_bet_2017_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2017_allsites.saf.idx bet_2019_allsites.saf.idx -P 1 > bet_2017_2019_allsites.sfs ; realSFS fst index bet_2017_allsites.saf.idx bet_2019_allsites.saf.idx -sfs bet_2017_2019_allsites.sfs -fstout bet_2017_2019_allsites

qsub sfs_bet_2017_2019

# Global Fst between populations 
realSFS fst stats bet_2017_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:99207]:0.390559 Fst.Weight:0.246598


# Writing down 2d-SFS priors - 2017 & 2020
>sfs_bet_2017_2020
nano sfs_bet_2017_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2017_allsites.saf.idx bet_2020_allsites.saf.idx -P 1 > bet_2017_2020_allsites.sfs ; realSFS fst index bet_2017_allsites.saf.idx bet_2020_allsites.saf.idx -sfs bet_2017_2020_allsites.sfs -fstout bet_2017_2020_allsites

qsub sfs_bet_2017_2020

# Global Fst between populations 
realSFS fst stats bet_2017_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:158601]:0.028936 Fst.Weight:0.084642


# Writing down 2d-SFS priors - 2017 & 2022
>sfs_bet_2017_2022
nano sfs_bet_2017_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2017_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2017_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2017_allsites.saf.idx bet_2022_allsites.saf.idx -P 1 > bet_2017_2022_allsites.sfs ; realSFS fst index bet_2017_allsites.saf.idx bet_2022_allsites.saf.idx -sfs bet_2017_2022_allsites.sfs -fstout bet_2017_2022_allsites

qsub sfs_bet_2017_2022

# Global Fst between populations 
realSFS fst stats bet_2017_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:192714]:0.298955 Fst.Weight:0.085998


# Writing down 2d-SFS priors - 2018 & 2019
>sfs_bet_2018_2019
nano sfs_bet_2018_2019

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018_2019 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018_2019.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2018_allsites.saf.idx bet_2019_allsites.saf.idx -P 1 > bet_2018_2019_allsites.sfs ; realSFS fst index bet_2018_allsites.saf.idx bet_2019_allsites.saf.idx -sfs bet_2018_2019_allsites.sfs -fstout bet_2018_2019_allsites

qsub sfs_bet_2018_2019

# Global Fst between populations 
realSFS fst stats bet_2018_2019_allsites.fst.idx

# Output:
FST.Unweight[nObs:98142]:0.530621 Fst.Weight:0.177405


# Writing down 2d-SFS priors - 2018 & 2020
>sfs_bet_2018_2020
nano sfs_bet_2018_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2018_allsites.saf.idx bet_2020_allsites.saf.idx -P 1 > bet_2018_2020_allsites.sfs ; realSFS fst index bet_2018_allsites.saf.idx bet_2020_allsites.saf.idx -sfs bet_2018_2020_allsites.sfs -fstout bet_2018_2020_allsites

qsub sfs_bet_2018_2020

# Global Fst between populations 
realSFS fst stats bet_2018_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:157872]:0.009996 Fst.Weight:0.029215


# Writing down 2d-SFS priors - 2018 & 2022
>sfs_bet_2018_2022
nano sfs_bet_2018_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2018_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2018_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2018_allsites.saf.idx bet_2022_allsites.saf.idx -P 1 > bet_2018_2022_allsites.sfs ; realSFS fst index bet_2018_allsites.saf.idx bet_2022_allsites.saf.idx -sfs bet_2018_2022_allsites.sfs -fstout bet_2018_2022_allsites

qsub sfs_bet_2018_2022

# Global Fst between populations 
realSFS fst stats bet_2018_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:187940]:0.017806 Fst.Weight:0.030153


# Writing down 2d-SFS priors - 2019 & 2020
>sfs_bet_2019_2020
nano sfs_bet_2019_2020

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2019_2020 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2019_2020.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2019_allsites.saf.idx bet_2020_allsites.saf.idx -P 1 > bet_2019_2020_allsites.sfs ; realSFS fst index bet_2019_allsites.saf.idx bet_2020_allsites.saf.idx -sfs bet_2019_2020_allsites.sfs -fstout bet_2019_2020_allsites

qsub sfs_bet_2019_2020

# Global Fst between populations 
realSFS fst stats bet_2019_2020_allsites.fst.idx

# Output:
FST.Unweight[nObs:88003]:0.349794 Fst.Weight:0.166096


# Writing down 2d-SFS priors - 2019 & 2022
>sfs_bet_2019_2022
nano sfs_bet_2019_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2019_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2019_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2019_allsites.saf.idx bet_2022_allsites.saf.idx -P 1 > bet_2019_2022_allsites.sfs ; realSFS fst index bet_2019_allsites.saf.idx bet_2022_allsites.saf.idx -sfs bet_2019_2022_allsites.sfs -fstout bet_2019_2022_allsites

qsub sfs_bet_2019_2022

# Global Fst between populations 
realSFS fst stats bet_2019_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:99359]:0.708405 Fst.Weight:0.185706


# Writing down 2d-SFS priors - 2020 & 2022
>sfs_bet_2020_2022
nano sfs_bet_2020_2022

#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N sfs_bet_2020_2022 # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M jaskielj@bu.edu #your email
#$ -m bea
#$ -pe omp 1
#$ -j y # Join standard output and error to a single file
#$ -o sfs_bet_2020_2022.qlog # Name the file where to redirect standard output and error

module load angsd
realSFS bet_2020_allsites.saf.idx bet_2022_allsites.saf.idx -P 1 > bet_2020_2022_allsites.sfs ; realSFS fst index bet_2020_allsites.saf.idx bet_2022_allsites.saf.idx -sfs bet_2020_2022_allsites.sfs -fstout bet_2020_2022_allsites

qsub sfs_bet_2020_2022

# Global Fst between populations 
realSFS fst stats bet_2020_2022_allsites.fst.idx

# Output:
FST.Unweight[nObs:158859]:-0.001093 Fst.Weight:0.036291


#------------------------------------------ Tajima's D ------------------------------------------
##ANGSD EXAMPLE
#(from https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests)
module load angsd
module list
#angsd/0.935

#Step 1: Finding a 'global estimate' of the SFS
#First estimate the site allele frequency likelihood
./angsd -bam bam.filelist -doSaf 1 -anc chimpHg19.fa -GL 1 -P 24 -out out

#Obtain the maximum likelihood estimate of the SFS using the realSFS program found in the misc subfolder
./misc/realSFS out.saf.idx -P 24 > out.sfs

#You did this above when you were calculating SAF and SFS for each pop group (use filters that don't alter allele frequency for this)
#Trying folded and unfolded to see difference (spoiler, extremely little difference)
realSFS skj_dem.saf.idx -fold 1 > skj_dem_folded.sfs
realSFS yft_dem.saf.idx -fold 1 > yft_dem_folded.sfs
realSFS bet_dem.saf.idx -fold 1 > bet_dem_folded.sfs

realSFS skj_dem.saf.idx > skj_dem.sfs
realSFS yft_dem.saf.idx > yft_dem.sfs
realSFS bet_dem.saf.idx > bet_dem.sfs

#Step 2: Calculate the thetas for each site
realSFS saf2theta out.saf.idx -sfs out.sfs -outname out

#For yours:
#Folded
realSFS saf2theta skj_dem.saf.idx -sfs skj_dem_folded.sfs -outname skj_theta_noLD_fold -fold 1
realSFS saf2theta yft_dem.saf.idx -sfs yft_dem_folded.sfs -outname yft_theta_noLD_fold -fold 1
realSFS saf2theta bet_dem.saf.idx -sfs bet_dem_folded.sfs -outname bet_theta_noLD_fold -fold 1
#Unfolded
realSFS saf2theta skj_dem.saf.idx -sfs skj_dem.sfs -outname skj_theta_noLD
realSFS saf2theta yft_dem.saf.idx -sfs yft_dem.sfs -outname yft_theta_noLD
realSFS saf2theta bet_dem.saf.idx -sfs bet_dem.sfs -outname bet_theta_noLD

#The output from the above command are two files out.thetas.gz and out.thetas.idx. A formal description of these files can be found in the doc/formats.pdf in the angsd package. It is possible to extract the logscale persite thetas using the ./thetaStat print program.
thetaStat print out.thetas.idx 2>/dev/null |head

#Mine:
#Folded
thetaStat print skj_theta_noLD_fold.thetas.idx 2>/dev/null |head
thetaStat print yft_theta_noLD_fold.thetas.idx 2>/dev/null |head 
thetaStat print bet_theta_noLD_fold.thetas.idx 2>/dev/null |head
#Unfolded
thetaStat print skj_theta_noLD.thetas.idx 2>/dev/null |head
thetaStat print yft_theta_noLD.thetas.idx 2>/dev/null |head 
thetaStat print bet_theta_noLD.thetas.idx 2>/dev/null |head

#Step 3a: Estimate Tajimas D and other statistics
#calculate Tajimas D
./misc/thetaStat do_stat out.thetas.idx

#For mine:
#Folded
#skj 
thetaStat do_stat skj_theta_noLD_fold.thetas.idx
cat skj_theta_noLD_fold.thetas.idx.pestPG

#yft
thetaStat do_stat yft_theta_noLD_fold.thetas.idx
cat yft_theta_noLD_fold.thetas.idx.pestPG

#bet
thetaStat do_stat bet_theta_noLD_fold.thetas.idx
cat bet_theta_noLD_fold.thetas.idx.pestPG

#Unfolded
#skj 
thetaStat do_stat skj_theta_noLD.thetas.idx
cat skj_theta_noLD.thetas.idx.pestPG

#yft
thetaStat do_stat yft_theta_noLD.thetas.idx
cat yft_theta_noLD.thetas.idx.pestPG

#bet
thetaStat do_stat bet_theta_noLD.thetas.idx
cat bet_theta_noLD.thetas.idx.pestPG

#Step 3b: Sliding Window example
#We can easily do a sliding window analysis by adding -win/-step arguments to the last command
#This will calculate the test statistic using a window size of 50kb and a step size of 10kb.
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz

#For mine:
#Folded
#skj ref
thetaStat do_stat skj_theta_noLD_fold.thetas.idx -win 50000 -step 10000 -outnames skj_noLD_fold_50kbwin10kbstep.thetasWindow.gz
#output in skj_noLD_fold_50kbwin10kbstep.thetasWindow.gz.pestPG
#yft
thetaStat do_stat yft_theta_noLD_fold.thetas.idx -win 50000 -step 10000 -outnames yft_noLD_fold_50kbwin10kbstep.thetasWindow.gz
#output in yft_noLD_fold_50kbwin10kbstep.thetasWindow.gz.pestPG
#bet
thetaStat do_stat bet_theta_noLD_fold.thetas.idx -win 50000 -step 10000 -outnames bet_noLD_fold_50kbwin10kbstep.thetasWindow.gz
##output in bet_noLD_fold_50kbwin10kbstep.thetasWindow.gz.pestPG

#Unfolded
#skj ref
thetaStat do_stat skj_theta_noLD.thetas.idx -win 50000 -step 10000 -outnames skj_noLD_50kbwin10kbstep.thetasWindow.gz
#output in skj_noLD_50kbwin10kbstep.thetasWindow.gz.pestPG
#yft
thetaStat do_stat yft_theta_noLD.thetas.idx -win 50000 -step 10000 -outnames yft_noLD_50kbwin10kbstep.thetasWindow.gz
#output in yft_noLD_50kbwin10kbstep.thetasWindow.gz.pestPG
#bet
thetaStat do_stat bet_theta_noLD.thetas.idx -win 50000 -step 10000 -outnames bet_noLD_50kbwin10kbstep.thetasWindow.gz
##output in bet_noLD_50kbwin10kbstep.thetasWindow.gz.pestPG

#- Output in the ./thetaStat print thetas.idx are the log scaled per site estimates of the thetas
#- Output in the pestPG file are the sum of the per site estimates for a region


#The .pestPG file is a 14 column file (tab seperated). The first column contains information about the region. The second and third column is the reference name and the center of the window.
#We then have 5 different estimators of theta, these are: Watterson, pairwise, FuLi, fayH, L. And we have 5 different neutrality test statistics: Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E. The final column is the effetive number of sites with data in the window.
#Format is:
#(indexStart,indexStop)(posStart,posStop)(regStat,regStop) chrname wincenter tW tP tF tH tL tajD fulif fuliD fayH zengsE numSites
#Most likely you are just interest in the wincenter (column 3) and the column 9 which is the Tajima's D statistic.
#The first 3 columns relates to the region. The next 5 columns are 5 different estimators of theta, and the next 5 columns are neutrality test statistics. The final column is the number of sites with data in the region.

#See R doc for plotting
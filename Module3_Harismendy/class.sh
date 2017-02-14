#!/bin/bash

################################################################
#
#  "Looking at raw data: reads, alignment and variants"
#	MED263 Harismendy 
#	W2017
#
################################################################


###########
#
# PreRequisite
#
###########

#create a local variable for the path to the class dir
classdir=/projects/ps-yeolab/biom262_2017/harismendy/MED263_harismendy

#in your home directory create a workingdir for the class
mkdir workingdir 
cd workingdir 	
	

###########
#
# Using screen
#
###########

		#Use screen
		screen –S name to create
		
		#^A ^D to detach
		screen –x [PID].name #to attach
		exit #to kill
		
###########
#
# UNIX/AWK to work with intervals
#
###########
		
#How many lines in the CGC.bed file? 
	wc –l $classdir/resources/CGC.exons.bed
#How many CGC exons on chromosome 1 ?
	awk ‘$1==“chr1”’ $classdir/resources/CGC.exons.bed | wc -l
#how many CGC genes on chromsome 1? 
	cut –f 4 $classdir/resources/CGC.exons.bed | sort | uniq | wc –l
#how many exons per CGC gene ? 
	cut –f 4 $classdir/resources/CGC.exons.bed | sort | uniq –c | sort –r | more
#what is the longest gene in the list ? 
	awk ‘$5=$3-$2’ $classdir/resources/CGC.exons.bed | sort –nrk 5 | head
#what is the sum of the length of all the exons
	awk ‘{sum=sum+$3-$2} END {print sum/NR}’ $classdir/resources/CGC.exons.bed

		
###########
#
# FASTQC and MULTIQC
#
###########

#Inspect the fastq file. How long are the reads? 
	zcat $classdir/materials/SRR866442_1.fastq.gz | more
#Create an output directory
	mkdir $workingdir/fastqc
#Run fastqc on all files
	for file in $classdir/materials/*SRR*fastq.gz; do fastqc –o $workingdir/fastqc $file & done
#Running MultiQC
	cd $workingdir/fastqc
	multiqc . # the dot indcates the current directory
	cd ..
	tar -cvf fastqc_results.tar fastqc #creates an archive of the results
#On your Desktop/Laptop
	#copy the archive to your latop
	#open the multiqc_report.html with your browser

###########
#
# Read Alignment (targeted amplicons)
#
###########
#create a bam folder
	mkdir $workingdir/bam
	cd $workingdir/bam

#Start a Screen session
	screen –S alignment

#Alignment + convert to sorted bam
	bwa mem $classdir/resources/chr1.fa.gz $classdir/materials/SRR866442_1.fastq.gz $classdir/materials/SRR866442_2.fastq.gz | samtools view –buSh - > SRR866442.bam

#Sort and index the bam file
	samtools sort –m 2G SRR866442.bam > SRR866442.sorted.bam
	samtools index SRR866442.sorted.bam
	
#What fraction of reads are properly aligned? 
	samtools flagstat SRR866442.sorted.bam > SRR866442.flagstat.txt


###########
#
# Statistics and Slices (whole exome)
#
###########

#What fraction of reads are properly aligned? 
	samtools flagstat $classdir/materials/CPTRES7.chr21.bam > CPTRES7.flagstat.txt

#How many “gapped reads” ?
	samtools view $classdir/materials/CPTRES7.chr21.bam | awk ‘$6~/[ID]/’ | wc -l

#how many reads over the CGC exons
	samtools view –L $classdir/resources/CGC.exons.bed $classdir/materials/CPTRES7.chr21.bam | wc -l

###########
#
# Visualize Alignments in IGV
#
###########

# 	1 copy the bam and bai file to your laptop
#	2. open with IGV
#


###########
#
# Calculating Coverage Depth
#
###########

#what fraction of the CGC chr1 are covered by more than 20 reads
	grep '^chr21' $classdir/resources/CGC.exons.bed | samtools depth -b - $classdir/materials/CPTRES7.chr21.bam | awk '$3>20' | wc -l  
	grep '^chr21' $classdir/resources/CGC.exons.bed | awk '{sum+=$3-$2} END {print sum}' #total number of CGC exons bp on chr21

#how many RUNX1 base pairs are covered at 20x or greater? 
	bedtools coverage -a $classdir/resourcesCGC.exons.bed -b $classdir/material/CPTRES7.chr21.bam -hist | grep 'RUNX1' | awk '$5>20' | awk '{sum+=$6} END {print sum}' #solution #1
	grep 'RUNX1' $classdir/resources/CGC.exons.bed | samtools depth -b - $classdir/materials/CPTRES7.chr21.bam | awk '$3>20' | wc -l #solution #2


###########
#
# Filtering Variants
#
###########
#create variant directory
mkdir $workingdir/variants

# index the vcf file (will create in the classdir folder. only one needed !)
 	tabix -p vcf $classdir/materials/CPTRES1vs15.vcf.gz
 	
# how many variants pass the filters?
	zgrep -v '^##' $classdir/materials/CPTRES1vs15.vcf.gz | awk '$7=="PASS"' | wc -l
 	
# Filter the variants 
 	bcftools filter -i 'FILTER=="PASS"' $classdir/materials/CPTRES1vs15.vcf.gz > CPTRES1vs15.PASS.vcf.gz
 
# What is the transition to transversion ratio? 
	bcftools stats CPTRES1vs15.PASS.vcf.gz
	
# how many somatic variants ?
	bcftools filter -i 'INFO/SS==2' CPTRES1vs15.PASS.vcf.gz | wc -l

# export the soamtic variant to tsv
	bcftools filter -i 'INFO/SS==2' CPTRES1vs15.PASS.vcf.gz | bcftools view -O z > CPTRES1vs15.PASS.SOM.vcf.gz
	vcf2tsv -g CPTRES1vs15.PASS.SOM.vcf.gz > CPTRES1vs15.PASS.SOM.tsv

# How many nonsense soamtic variants

	#annotate with annovar 
	$classdir/resources/annovar/table_annovar.pl --vcfinput --nastring . --protocol refGene --operation g --buildver hg19 --outfile CPTRESann CPTRES1vs15.PASS.SOM.vcf.gz $classdir/resources/annovar/humandb/
	# count the 2nd column of the variant_function file
	cut -f 2 CPTRESann.refGene.exonic_variant_function | sort | uniq -c

###########
#
# Using R for data wrangling and plotting 
#
###########

# Import variant TSV into R


	# aggregate, annotate
	#going over dplyr, reshape, aggregate

	# plot
	# going over ggplot






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
#


classdir=/path/to/class/directory # the class directory cotnains both resoruces and material directory


mkdir harismendy_workingdir # the following command assume workingdir is inside the MED263_harismendy directory
cd harismendy_workingdir 	
	

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
#what is the gene with most exons in the list ? 
	awk ‘$5=$3-$2’ $classdir/resources/CGC.exons.bed | sort –nrk 5 | head
#what is the sum of the length of all the exons
	awk ‘{sum=sum+$3-$2} END {print sum/NR}’ $classdir/resources/CGC.exons.bed

		
###########
#
# FASTQC and MULTIQC
#
###########
# create an fastqc subdirectory

mkdir fastqc
cd fastqc

#Inspect the fastq file. How long are the reads? 
	zcat $classdir/materials/SRR866442_1.fastq.gz | more
#Create an output directory
	mkdir fastqc
#Run fastqc on all files
	for file in $classdir/materials/*SRR*fastq.gz; do fastqc –o fastqc $file & done
#Running MultiQC
	cd fastqc
	multiqc . # the dot indcates the current directory
	cd ..
	tar -cvf fastqc_results.tar fastqc #creates an archive of the results
#On your Desktop/Laptop
	scp -P 9221 username@server:/home/username/harismendy_workingdir/fastqc_results.tar . # copy the results
	#open the multiqc_report.html with your browser

###########
#
# Read Alignment (targeted amplicons)
#
###########


# create an alignment subdirectory
cd harismendy_workingdir
mkdir alignment 
cd alignment

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
# Statistics and Slices (whole exome, chr21 slice)
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
	bedtools coverage -a $classdir/materials/CPTRES7.chr21.bam -b $classdir/resources/CGC.exons.bed  -hist | grep 'RUNX1' | awk '$5>20' | awk '{sum+=$6} END {print sum}' #solution #1
	grep 'RUNX1' $classdir/resources/CGC.exons.bed | samtools depth -b - $classdir/materials/CPTRES7.chr21.bam | awk '$3>20' | wc -l #solution #2


###########
#
# Filtering Variants
#
###########

# create an variants subdirectory
cd harismendy_workingdir
mkdir variants 
cd variants

# index the vcf file
 	tabix -p vcf $classdir/materials/CPTRES1vs15.vcf.gz
 	
# how many variants pass the filters?
	zgrep -v '^##' $classdir/materials/CPTRES1vs15.vcf.gz | awk '$7=="PASS"' | wc -l
 	
# Filter the variants 
 	bcftools filter -i 'FILTER=="PASS"' $classdir/materials/CPTRES1vs15.vcf.gz | bcftools view -O z > CPTRES1vs15.PASS.vcf.gz
 
# What is the transition to transversion ratio? 
	bcftools stats CPTRES1vs15.PASS.vcf.gz
	
# how many somatic variants ?
	bcftools filter -i 'INFO/SS==2' CPTRES1vs15.PASS.vcf.gz > CPTRES1vs15.PASS.SOM.vcf 
	grep -v '^#' CPTRES1vs15.PASS.SOM.vcf | wc -l

# how many LOH variants per CGC exons? 
	

# export the somatic variants (mutation and LOH) to tsv
	bcftools filter -i 'INFO/SS==2 | INFO/SS==3' CPTRES1vs15.PASS.SOM.vcf.gz | bcftools view -O z > CPTRES1vs15.PASS.SOM.vcf.gz
	vcf2tsv -g CPTRES1vs15.PASS.vcf.gz > CPTRES1vs15.PASS.SOM.tsv

# How many nonsense somatic variants?

	#annotate with annovar 
	$classdir/resources/annovar/table_annovar.pl --vcfinput --nastring . --protocol refGene --operation g --buildver hg19 --outfile CPTRESann CPTRES1vs15.PASS.vcf.gz $classdir/resources/annovar/humandb/
	# count the 2nd column of the variant_function file
	cut -f 2 CPTRESann.refGene.exonic_variant_function | sort | uniq -c
	
	

###########
#
# Using RStudio and R for data wrangling and plotting 
#
###########

#load the dplyr and ggplot2, reshape and ggplot library

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

# Import variant TSV into R

	variants<-read.delim("CPTRES1vs15.PASS.all.tsv",header=TRUE)


# number of unique variants

	nrow(unique(select(variants,CHROM,POS,REF,ALT)))
	variants %>% select(CHROM,POS,REF,ALT) %>% unique() %>% tally()

# number of variants in each SS category

	variants %>% select(CHROM,POS,REF,ALT,SS) %>% unique() %>% group_by(SS) %>% tally()

# number of variants in SS x VT categories

	variants %>% select(CHROM,POS,REF,ALT,SS,VT) %>% unique() %>% group_by(SS,VT) %>% tally()

# distribution of GT value in the control sample

	variants %>% filter(SAMPLE=="control") %>% select(CHROM,POS,REF,ALT,GT) %>% unique() %>% group_by(GT) %>% tally()


#calculate the allele fraction from the AD field
	variants$AD<-as.character(variants$AD)
	variants<-variants %>% separate(AD,c("AD1","AD2"),sep=",")
	
	variants$AD1<-as.numeric(variants$AD1)
	variants$AD2<-as.numeric(variants$AD2)
	variants<-mutate(variants,AFrac=(AD2/(AD1+AD2)))

#Density plot of AFrac as a function of SS in the treated
	ggplot(filter(variants,SAMPLE=="treated"),aes(AFrac,fill=factor(SS)))+geom_density(alpha=0.5)

#boxplot of DP vs SS
	ggplot(variants,aes(factor(SS),log10(DP)))+geom_boxplot() #overall
	ggplot(variants,aes(factor(SS),log10(AD1+AD2),fill=SAMPLE))+geom_boxplot() #per sample
	
#Scatter plots of AFrac in control vs Treated

	#Cast the AFrac values into two column control treated
	FracCast<-dcast(CHROM+POS+ID+REF+ALT+VT+SS~SAMPLE,data=variants,value.var="AFrac")

	#Scatter plot  Add VT, CHROM, SS information through colors or shapes. 

	ggplot(FracCast,aes(control,treated,col=factor(SS)))+geom_point(alpha=0.5)
	
	


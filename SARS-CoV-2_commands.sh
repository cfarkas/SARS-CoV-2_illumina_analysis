#!/bin/bash

{

#########################
#########################
### Illumina Analysis ###
#########################
#########################

#############
# Downloads #
#############

### March 20, 2020 Paired end and Single End


fastq-dump -Z SRR10903401 > SRR10903401.fq
fastq-dump -Z SRR10903402 > SRR10903402.fq
fastq-dump -Z SRR10971381 > SRR10971381.fq
fastq-dump -Z SRR11059940 > SRR11059940.fq
fastq-dump -Z SRR11059941 > SRR11059941.fq
fastq-dump -Z SRR11059942 > SRR11059942.fq
fastq-dump -Z SRR11059943 > SRR11059943.fq
fastq-dump -Z SRR11059944 > SRR11059944.fq
fastq-dump -Z SRR11059945 > SRR11059945.fq
fastq-dump -Z SRR11059946 > SRR11059946.fq
fastq-dump -Z SRR11059947 > SRR11059947.fq
fastq-dump -Z SRR11140744 > SRR11140744.fq
fastq-dump -Z SRR11140746 > SRR11140746.fq
fastq-dump -Z SRR11140748 > SRR11140748.fq
fastq-dump -Z SRR11140750 > SRR11140750.fq
fastq-dump -Z SRR11177792 > SRR11177792.fq
fastq-dump -Z SRR11241254 > SRR11241254.fq
fastq-dump -Z SRR11241255 > SRR11241255.fq
fastq-dump -Z SRR11247075 > SRR11247075.fq
fastq-dump -Z SRR11247076 > SRR11247076.fq
fastq-dump -Z SRR11247077 > SRR11247077.fq
fastq-dump -Z SRR11247078 > SRR11247078.fq
fastq-dump -Z SRR11278090 > SRR11278090.fq
fastq-dump -Z SRR11278091 > SRR11278091.fq
fastq-dump -Z SRR11278092 > SRR11278092.fq
fastq-dump -Z SRR11278164 > SRR11278164.fq
fastq-dump -Z SRR11278165 > SRR11278165.fq
fastq-dump -Z SRR11278166 > SRR11278166.fq
fastq-dump -Z SRR11278167 > SRR11278167.fq
fastq-dump -Z SRR11278168 > SRR11278168.fq
fastq-dump -Z SRR11314339 > SRR11314339.fq

### March 25, 2020 Paired end and Single End

fastq-dump -Z SRR11397714 > SRR11397714.fq
fastq-dump -Z SRR11397715 > SRR11397715.fq
fastq-dump -Z SRR11397716 > SRR11397716.fq
fastq-dump -Z SRR11397717 > SRR11397717.fq
fastq-dump -Z SRR11397718 > SRR11397718.fq
fastq-dump -Z SRR11397719 > SRR11397719.fq
fastq-dump -Z SRR11397720 > SRR11397720.fq
fastq-dump -Z SRR11397721 > SRR11397721.fq
fastq-dump -Z SRR11397728 > SRR11397728.fq
fastq-dump -Z SRR11397729 > SRR11397729.fq
fastq-dump -Z SRR11397730 > SRR11397730.fq
fastq-dump -Z SRR11393704 > SRR11393704.fq


### Bowtie commands
bowtie2-build covid19-refseq.fasta covid19-refseq
samtools faidx covid19-refseq.fasta

# Alignment
a= ls -1 *.fq
for a in *.fq; do bowtie2 -p 20 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq ${a} > ${a}.sam
done

# BAM file generation
b= ls -1 *.sam
for b in *.sam; do samtools sort ${b} > ${b}.sorted.bam -@ 20
done

# Remove duplicates
c= ls -1 *.fq.sam.sorted.bam
for c in *.fq.sam.sorted.bam; do samtools rmdup ${c} ${c}.rmdup
done

# Index bam files
d= ls -1 *.fq.sam.sorted.bam
for d in *.fq.sam.sorted.bam; do samtools index ${d} -@ 20
done

# Variant Calling
e= ls -1 *.fq.sam.sorted.bam.rmdup
for e in *.fq.sam.sorted.bam.rmdup; do bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou ${e}| bcftools call -mv -Ov -o ${e}.vcf
done
rm *.sam

# VCF filtering
f= ls -1 *.fq.sam.sorted.bam.rmdup.vcf
for f in *.fq.sam.sorted.bam.rmdup.vcf; do bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${f} > ${f}.filtered
done

mkdir filtered_vcfs
cp *.filtered ./filtered_vcfs/

# Extracting DP4 field to extract allele frequency of variants

cd filtered_vcfs
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR10903401.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR10903401.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR10903402.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR10903402.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR10971381.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR10971381.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059940.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059940.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059941.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059941.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059942.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059942.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059943.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059943.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059944.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059944.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059945.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059945.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059946.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059946.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11059947.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11059947.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11140744.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11140744.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11140746.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11140746.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11140748.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11140748.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11140750.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11140750.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11177792.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11177792.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11241254.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11241254.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11241255.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11241255.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11247075.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11247075.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11247076.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11247076.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11247077.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11247077.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11247078.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11247078.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278090.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278090.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278091.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278091.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278092.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278092.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278164.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278164.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278165.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278165.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278166.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278166.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278167.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278167.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11278168.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11278168.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11314339.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11314339.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397714.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397714.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397715.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397715.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397716.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397716.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397717.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397717.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397718.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397718.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397719.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397719.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397720.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397720.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397721.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397721.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397728.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397728.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397729.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397729.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11397730.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11397730.DP4
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' SRR11393704.fq.sam.sorted.bam.rmdup.vcf.filtered > SRR11393704.DP4

cd ..

#####################
# Coverage Analysis #
#####################

a= ls -1 *.fq.sam.sorted.bam
for a in *.fq.sam.sorted.bam; do samtools depth ${a} |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
done

############################################################################################################
############################################################################################################

#######################
#######################
### Snippy Analysis ###
#######################
#######################

mkdir Snippy_results
cd Snippy_results

# China

snippy --cpus 20 --outdir ./SRR11059940 --ref SARS-CoV-19.gb --se SRR11059940.fq
snippy --cpus 20 --outdir ./SRR11059944 --ref SARS-CoV-19.gb --se SRR11059944.fq
snippy --cpus 20 --outdir ./SRR11059945 --ref SARS-CoV-19.gb --se SRR11059945.fq
snippy --cpus 20 --outdir ./SRR11059946 --ref SARS-CoV-19.gb --se SRR11059946.fq
snippy --cpus 20 --outdir ./SRR11059947 --ref SARS-CoV-19.gb --se SRR11059947.fq

# USA
snippy --cpus 20 --outdir ./SRR11241254 --ref SARS-CoV-19.gb --se SRR11241254.fq
snippy --cpus 20 --outdir ./SRR11241255 --ref SARS-CoV-19.gb --se SRR11241255.fq
snippy --cpus 20 --outdir ./SRR11247075 --ref SARS-CoV-19.gb --se SRR11247075.fq
snippy --cpus 20 --outdir ./SRR11247076 --ref SARS-CoV-19.gb --se SRR11247076.fq
snippy --cpus 20 --outdir ./SRR11247077 --ref SARS-CoV-19.gb --se SRR11247077.fq
snippy --cpus 20 --outdir ./SRR11247078 --ref SARS-CoV-19.gb --se SRR11247078.fq
snippy --cpus 20 --outdir ./SRR11278090 --ref SARS-CoV-19.gb --se SRR11278090.fq
snippy --cpus 20 --outdir ./SRR11278091 --ref SARS-CoV-19.gb --se SRR11278091.fq
snippy --cpus 20 --outdir ./SRR11278092 --ref SARS-CoV-19.gb --se SRR11278092.fq
snippy --cpus 20 --outdir ./SRR11278164 --ref SARS-CoV-19.gb --se SRR11278164.fq
snippy --cpus 20 --outdir ./SRR11278165 --ref SARS-CoV-19.gb --se SRR11278165.fq
snippy --cpus 20 --outdir ./SRR11278166 --ref SARS-CoV-19.gb --se SRR11278166.fq
snippy --cpus 20 --outdir ./SRR11278167 --ref SARS-CoV-19.gb --se SRR11278167.fq
snippy --cpus 20 --outdir ./SRR11278168 --ref SARS-CoV-19.gb --se SRR11278168.fq

# Australia
snippy --cpus 20 --outdir ./SRR11397714 --ref SARS-CoV-19.gb --se SRR11397714.fq
snippy --cpus 20 --outdir ./SRR11397715 --ref SARS-CoV-19.gb --se SRR11397715.fq
snippy --cpus 20 --outdir ./SRR11397716 --ref SARS-CoV-19.gb --se SRR11397716.fq
snippy --cpus 20 --outdir ./SRR11397717 --ref SARS-CoV-19.gb --se SRR11397717.fq
snippy --cpus 20 --outdir ./SRR11397718 --ref SARS-CoV-19.gb --se SRR11397718.fq
snippy --cpus 20 --outdir ./SRR11397719 --ref SARS-CoV-19.gb --se SRR11397719.fq
snippy --cpus 20 --outdir ./SRR11397720 --ref SARS-CoV-19.gb --se SRR11397720.fq
snippy --cpus 20 --outdir ./SRR11397721 --ref SARS-CoV-19.gb --se SRR11397721.fq
snippy --cpus 20 --outdir ./SRR11397728 --ref SARS-CoV-19.gb --se SRR11397728.fq
snippy --cpus 20 --outdir ./SRR11397729 --ref SARS-CoV-19.gb --se SRR11397729.fq
snippy --cpus 20 --outdir ./SRR11397730 --ref SARS-CoV-19.gb --se SRR11397730.fq

cd ..

################################
################################
### Alignment of CDC primers ###
################################
################################


bowtie2 -p 20 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f CDC_primers.fasta > CDC_primers.sam
samtools view -bS CDC_primers.sam > CDC_primers.bam
samtools sort -o CDC_primers.sorted.bam CDC_primers.bam
bedtools bamtobed -i CDC_primers.sorted.bam > CDC_primers.bed

####################################################################################################
####################################################################################################
### Alignment of Hong Kong University, Pasteur Institute Primers and Korean Primers (PMC7045880) ###
####################################################################################################
####################################################################################################

bowtie2 -p 20 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f HK_Pasteur_Korea.fasta > HK_Pasteur_Korea.sam
samtools view -bS HK_Pasteur_Korea.sam > HK_Pasteur_Korea.bam
samtools sort -o HK_Pasteur_Korea.sorted.bam HK_Pasteur_Korea.bam
bedtools bamtobed -i HK_Pasteur_Korea.sorted.bam > HK_Pasteur_Korea.bed

###########################################################################################
###########################################################################################
### Alignment of Primer-blast primers: https://www.ncbi.nlm.nih.gov/tools/primer-blast/ ###
###########################################################################################
###########################################################################################

bowtie2 -p 20 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f primer-blast-50-170-bp.fasta > primer-blast-50-170-bp.sam
samtools view -bS primer-blast-50-170-bp.sam > primer-blast-50-170-bp.bam
samtools sort -o primer-blast-50-170-bp.sorted.bam primer-blast-50-170-bp.bam
bedtools bamtobed -i primer-blast-50-170-bp.sorted.bam > primer-blast-50-170-bp.bed


######################################################
######################################################
### Merging all VCFs and intersecting with primers ###
######################################################
######################################################

a= ls -1 *.fq.sam.sorted.bam.rmdup.vcf.filtered
for a in *.fq.sam.sorted.bam.rmdup.vcf.filtered; do bgzip ${a}
done

b= ls -1 *.fq.sam.sorted.bam.rmdup.vcf.filtered.gz
for b in *.fq.sam.sorted.bam.rmdup.vcf.filtered.gz; do tabix -p vcf ${b}
done

vcf-merge $(ls -1 *.fq.sam.sorted.bam.rmdup.vcf.filtered.gz | perl -pe 's/\n/ /g') > merge.vcf

### Intersecting primers

bedtools intersect -a CDC_primers.bed -b merge.vcf > CDC_primers.intersection
bedtools intersect -a HK_Pasteur_Korea.bed -b merge.vcf > HK_Pasteur_Korea.intersection
bedtools intersect -a primer-blast-50-170-bp.bed -b merge.vcf > primer-blast-50-170-bp.intersection


######################################################
######################################################
### Alignment of Genbank Sequences, March 31, 2020 ###
######################################################
######################################################

bowtie2 -p 5 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f genbank_sequences_March_31_2020.fasta > genbank_sequences_alignment.sam
samtools view -bS genbank_sequences_alignment.sam > genbank_sequences_alignment.bam
samtools sort -o genbank_sequences_alignment.sorted.bam genbank_sequences_alignment.bam
bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou genbank_sequences_alignment.sorted.bam| bcftools call -mv -Ov -o genbank_sequences_alignment.sorted.bam.vcf
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' genbank_sequences_alignment.sorted.bam.vcf > genbank_sequences_alignment.sorted.bam.DP4


###################################################################
###################################################################
### Alignment of Genbank Sequences, March 31, 2020: USA samples ###
###################################################################
###################################################################

bowtie2 -p 5 -D 5 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f genbank_sequences_USA_March_31_2020.fasta > genbank_sequences_USA_alignment.sam
samtools view -bS genbank_sequences_USA_alignment.sam > genbank_sequences_USA_alignment.bam
samtools sort -o genbank_sequences_USA_alignment.sorted.bam genbank_sequences_USA_alignment.bam
bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou genbank_sequences_USA_alignment.sorted.bam| bcftools call -mv -Ov -o genbank_sequences_USA_alignment.sorted.bam.vcf
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' genbank_sequences_USA_alignment.sorted.bam.vcf > genbank_sequences_USA_alignment.sorted.bam.DP4


###############################################################
###############################################################
### Final Primer Intersections: Adding USA genbank datasets ###
###############################################################
###############################################################

cp genbank_sequences_USA_alignment.sorted.bam.vcf ./genbank_USA.vcf
bgzip genbank_USA.vcf
tabix -p vcf genbank_USA.vcf.gz
bgzip merge.vcf
tabix -p vcf merge.vcf.gz

vcf-merge merge.vcf.gz genbank_USA.vcf.gz > final_merge.vcf

bedtools intersect -a CDC_primers.bed -b final_merge.vcf > CDC_primers.intersection.final
bedtools intersect -a HK_Pasteur_Korea.bed -b final_merge.vcf > HK_Pasteur_Korea.intersection.final
bedtools intersect -a primer-blast-50-170-bp.bed -b final_merge.vcf > primer-blast-50-170-bp.intersection.final

###############################################################
#
} | tee logfile
#

#!/bin/bash

{

##############################################
##############################################
### Variant Calling from illumina datasets ###
##############################################
##############################################

echo "Downloading illumina datasets..."
echo ""

#############
# Downloads #
#############

### March 20, 2020 Paired end and Single End

echo "Downloading illumina SRA files from March 20 and 25, 2020"

prefetch -O ./ SRR10903401
prefetch -O ./ SRR10903402
prefetch -O ./ SRR10971381
prefetch -O ./ SRR11059940
prefetch -O ./ SRR11059941
prefetch -O ./ SRR11059942
prefetch -O ./ SRR11059943
prefetch -O ./ SRR11059944
prefetch -O ./ SRR11059945
prefetch -O ./ SRR11059946
prefetch -O ./ SRR11059947
prefetch -O ./ SRR11140744
prefetch -O ./ SRR11140746
prefetch -O ./ SRR11140748
prefetch -O ./ SRR11140750
prefetch -O ./ SRR11177792
prefetch -O ./ SRR11241254
prefetch -O ./ SRR11241255
prefetch -O ./ SRR11247075
prefetch -O ./ SRR11247076
prefetch -O ./ SRR11247077
prefetch -O ./ SRR11247078
prefetch -O ./ SRR11278090
prefetch -O ./ SRR11278091
prefetch -O ./ SRR11278092
prefetch -O ./ SRR11278164
prefetch -O ./ SRR11278165
prefetch -O ./ SRR11278166
prefetch -O ./ SRR11278167
prefetch -O ./ SRR11278168
prefetch -O ./ SRR11314339

### March 25, 2020 Paired end and Single End

prefetch -O ./ SRR11397714
prefetch -O ./ SRR11397715
prefetch -O ./ SRR11397716
prefetch -O ./ SRR11397717
prefetch -O ./ SRR11397718
prefetch -O ./ SRR11397719
prefetch -O ./ SRR11397720
prefetch -O ./ SRR11397721
prefetch -O ./ SRR11397728
prefetch -O ./ SRR11397729
prefetch -O ./ SRR11397730
#
prefetch -O ./ SRR11410528
prefetch -O ./ SRR11410529
prefetch -O ./ SRR11410532
prefetch -O ./ SRR11410533
prefetch -O ./ SRR11410536
prefetch -O ./ SRR11410538
prefetch -O ./ SRR11410540
prefetch -O ./ SRR11410541
prefetch -O ./ SRR11410542

echo "Done"
echo ""
echo "Converting SRA files to fastq.gz"
SRA= ls -1 *.sra
for SRA in *.sra; do fastq-dump --gzip ${SRA}
done

##################################################################################
# Trimming downloaded Illumina datasets with fastp, using 16 threads (-w option) #
##################################################################################

echo "Trimming downloaded Illumina datasets with fastp."
echo ""

a= ls -1 *.fastq.gz
for a in *.fastq.gz; do fastp -w 16 -i ${a} -o ${a}.fastp
done

###########################################################################################
# Aligning illumina datasets againts reference with minimap, using 20 threads (-t option) #
###########################################################################################

echo "Aligning illumina datasets againts reference with minimap, using 20 threads."
echo ""

b= ls -1 *.fastq.gz.fastp
for b in *.fastq.gz.fastp; do minimap2 -ax sr covid19-refseq.fasta ${b} > ${b}.sam -t 20
done

###################################################
# Sorting SAM files, using 20 threads (-@ option) #
###################################################

echo "Sorting SAM files, using 20 threads."
echo ""

c= ls -1 *.sam
for c in *.sam; do samtools sort ${c} > ${c}.sorted.bam -@ 20
done

#########################################
# Calculating coverage of aligned reads #
#########################################

echo "Calculating coverage of aligned reads."
echo ""

d= ls -1 *.sorted.bam
for d in *.sorted.bam; do samtools depth ${d} |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
done

######################
# Renaming BAM files #
######################

echo "Renaming files in bash"
for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastq.gz.fastp.sam.sorted//g')";  done

######################
# Indexing BAM files #
######################

echo "Indexing BAM files."
echo ""

f= ls -1 *.bam
for f in *.bam; do samtools index ${f}; done

###############################################################
### Performing Germline Variant Calling with strelka v2.9.2 ###
###############################################################

echo "Performing Variant Calling with strelka v2.9.2:"
echo ""
echo "for documentation, please see: https://github.com/Illumina/strelka"
echo ""
echo "downloading strelka binary from github repository"
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2

echo "Calling Founder Variants"
bam= ls -1 *.bam
for bam in *.bam; do 
./strelka-2.9.2.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam SRR10971381.bam \
    --bam ${bam} \
    --referenceFasta covid19-refseq.fasta \
    --runDir ${bam}.founder
# execution on a single local machine with 20 parallel jobs
${bam}.founder/runWorkflow.py -m local -j 20
cp ./${bam}.founder/results/variants/variants.vcf.gz ./strelka_germline_variants.vcf.gz
bgzip -d strelka_germline_variants.vcf.gz
grep "#" strelka_germline_variants.vcf > strelka_germline_variants_header.vcf
grep "PASS" strelka_germline_variants.vcf > strelka_germline_variants_PASS.vcf
grep -v "NoPassedVariantGTs" strelka_germline_variants_PASS.vcf > strelka_germline_variants_PASS2.vcf
rm strelka_germline_variants_PASS.vcf
cat strelka_germline_variants_header.vcf strelka_germline_variants_PASS2.vcf > strelka_germline_variants.filtered.vcf
rm strelka_germline_variants_header.vcf strelka_germline_variants_PASS2.vcf
mv strelka_germline_variants.filtered.vcf ./${bam}.founder.vcf
rm strelka_germline_variants.vcf
rm -r -f ${bam}.founder
done
rm SRR10971381.bam.founder.vcf
#########################
# Cleaning Up SAM files #
#########################

echo "Cleaning up strelka directories and intermediate files"
rm *.sam

###############################################
### Merging founder variants using jacquard ###
###############################################
echo "Merging variants using jacquard"
echo ""
echo "for information, please see: https://jacquard.readthedocs.io/en/v0.42/overview.html#why-would-i-use-jacquard"
mkdir to_translate
cp *.founder.vcf ./to_translate/
cd to_translate
jacquard translate --force ./ translated_vcfs
cd ..
cp ./to_translate/translated_vcfs/* ./
mkdir primary_vcfs
mv *founder.vcf ./primary_vcfs/
mkdir to_merge
mv *.translatedTags.vcf ./to_merge/
cd to_merge
jacquard merge ./ merged.vcf
cd ..
cp ./to_merge/merged.vcf ./
echo "All done. Merged vcf is called merged.vcf and is located in current directory"

########################################################
########################################################
### Obtaining BED files from Primer Sets alignments  ###
########################################################
########################################################

echo "Indexing SARS-CoV-2 reference genome"
bowtie2-build covid19-refseq.fasta covid19-refseq
echo "Done"

###################################################
# Obtaining BED files from CDC primers alignments #
###################################################

echo "Obtaining BED files from CDC primers alignments."
echo ""
bowtie2 -p 20 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f CDC_primers.fasta > CDC_primers.sam
samtools view -bS CDC_primers.sam > CDC_primers.bam
samtools sort -o CDC_primers.sorted.bam CDC_primers.bam
bedtools bamtobed -i CDC_primers.sorted.bam > CDC_primers.bed

#######################################################################################################################
# Obtaining BED files from Hong Kong University, Pasteur Institute Primers and Korean Primers Alignments (PMC7045880) #
#######################################################################################################################

echo "Obtaining BED files from Hong Kong University, Pasteur Institute Primers and Korean Primers Alignments"
echo ""

bowtie2 -p 20 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f HK_Pasteur_Korea.fasta > HK_Pasteur_Korea.sam
samtools view -bS HK_Pasteur_Korea.sam > HK_Pasteur_Korea.bam
samtools sort -o HK_Pasteur_Korea.sorted.bam HK_Pasteur_Korea.bam
bedtools bamtobed -i HK_Pasteur_Korea.sorted.bam > HK_Pasteur_Korea.bed

##############################################################################################################
# Obtaining BED files from Primer-blast primer alignments : https://www.ncbi.nlm.nih.gov/tools/primer-blast/ #
##############################################################################################################

echo "Obtaining BED files from Primer-blast primer alignments : https://www.ncbi.nlm.nih.gov/tools/primer-blast/"
echo ""

bowtie2 -p 20 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -x covid19-refseq -f primer-blast-50-170-bp.fasta > primer-blast-50-170-bp.sam
samtools view -bS primer-blast-50-170-bp.sam > primer-blast-50-170-bp.bam
samtools sort -o primer-blast-50-170-bp.sorted.bam primer-blast-50-170-bp.bam
bedtools bamtobed -i primer-blast-50-170-bp.sorted.bam > primer-blast-50-170-bp.bed

### Intersecting primers
vcf2bed --deletions < merged.vcf > merged_deletions.bed
vcf2bed --snvs < merged.vcf > merged_snvs.bed
bedops --everything merged_{deletions,snvs}.bed > merged.bed
bedtools intersect -a CDC_primers.bed -b merged.bed > CDC_primers.intersection
bedtools intersect -a HK_Pasteur_Korea.bed -b merged.bed > HK_Pasteur_Korea.intersection
bedtools intersect -a primer-blast-50-170-bp.bed -b merged.bed > primer-blast-50-170-bp.intersection


##############################################################################################################
##############################################################################################################
### Alignment of Genbank Sequences, joining variants with illumina datasets and final primer intersections ###
##############################################################################################################
##############################################################################################################


##################################################
# Alignment of Genbank Sequences, March 31, 2020 #
##################################################

echo "Alignment of Genbank Sequences, March-31-2020"
echo ""

minimap2 -ax asm5 covid19-refseq.fasta genbank_sequences_March_31_2020.fasta > genbank_sequences_alignment.sam
samtools view -bS genbank_sequences_alignment.sam > genbank_sequences_alignment.bam
samtools sort -o genbank_sequences_alignment.sorted.bam genbank_sequences_alignment.bam
bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou genbank_sequences_alignment.sorted.bam| bcftools call -mv -Ov -o genbank_sequences_alignment.sorted.bam.vcf
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' genbank_sequences_alignment.sorted.bam.vcf > genbank_sequences_alignment.sorted.bam.DP4


###############################################################
# Alignment of Genbank Sequences, March 31, 2020: USA samples #
###############################################################

echo "Alignment of Genbank Sequences, March-31-2020: USA samples"
echo ""

minimap2 -ax asm5 covid19-refseq.fasta genbank_sequences_USA_March_31_2020.fasta > genbank_sequences_USA_alignment.sam
samtools view -bS genbank_sequences_USA_alignment.sam > genbank_sequences_USA_alignment.bam
samtools sort -o genbank_sequences_USA_alignment.sorted.bam genbank_sequences_USA_alignment.bam
bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou genbank_sequences_USA_alignment.sorted.bam| bcftools call -mv -Ov -o genbank_sequences_USA_alignment.sorted.bam.vcf
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' genbank_sequences_USA_alignment.sorted.bam.vcf > genbank_sequences_USA_alignment.sorted.bam.DP4

################################################################
# Alignment of Genbank Sequences, April 22, 2020: Asia samples # 
################################################################

echo "Alignment of Genbank Sequences, April-22-2020: Asia samples"
echo ""
gunzip genbank_sequences_Asia_April_22_2020.fasta.gz
minimap2 -ax asm5 covid19-refseq.fasta genbank_sequences_Asia_April_22_2020.fasta > genbank_sequences_Asia_22_2020_alignment.sam
samtools view -bS genbank_sequences_Asia_22_2020_alignment.sam > genbank_sequences_Asia_22_2020_alignment.bam
samtools sort -o genbank_sequences_Asia_22_2020_alignment.sorted.bam genbank_sequences_Asia_22_2020_alignment.bam
bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou genbank_sequences_Asia_22_2020_alignment.sorted.bam| bcftools call -mv -Ov -o genbank_sequences_Asia_22_2020_alignment.sorted.bam.vcf
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' genbank_sequences_Asia_22_2020_alignment.sorted.bam.vcf > genbank_sequences_Asia_22_2020_alignment.sorted.bam.DP4

#########################################################################
# Alignment of Genbank Sequences, April 22, 2020: North America samples # 
#########################################################################

echo "Alignment of Genbank Sequences, April-22-2020: North America samples"
echo ""
gunzip genbank_sequences_North_America_April_22_2020.fasta.gz
minimap2 -ax asm5 covid19-refseq.fasta genbank_sequences_North_America_April_22_2020.fasta > genbank_sequences_North_America_22_2020_alignment.sam
samtools view -bS genbank_sequences_North_America_22_2020_alignment.sam > genbank_sequences_North_America_22_2020_alignment.bam
samtools sort -o genbank_sequences_North_America_22_2020_alignment.sorted.bam genbank_sequences_North_America_22_2020_alignment.bam
bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou genbank_sequences_North_America_22_2020_alignment.sorted.bam| bcftools call -mv -Ov -o genbank_sequences_North_America_22_2020_alignment.sorted.bam.vcf
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' genbank_sequences_North_America_22_2020_alignment.sorted.bam.vcf > genbank_sequences_North_America_22_2020_alignment.sorted.bam.DP4

##################################################################
# Alignment of Genbank Sequences, April 22, 2020: Europe samples # 
##################################################################

echo "Alignment of Genbank Sequences, April-22-2020: Europe samples"
echo ""
gunzip genbank_sequences_Europe_April_22_2020.fasta.gz
minimap2 -ax asm5 covid19-refseq.fasta genbank_sequences_Europe_April_22_2020.fasta > genbank_sequences_Europe_22_2020_alignment.sam
samtools view -bS genbank_sequences_Europe_22_2020_alignment.sam > genbank_sequences_Europe_22_2020_alignment.bam
samtools sort -o genbank_sequences_Europe_22_2020_alignment.sorted.bam genbank_sequences_Europe_22_2020_alignment.bam
bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 10 -Ou genbank_sequences_Europe_22_2020_alignment.sorted.bam| bcftools call -mv -Ov -o genbank_sequences_Europe_22_2020_alignment.sorted.bam.vcf
bcftools query -f'[%CHROM\t%POS\t%DP4\n]' genbank_sequences_Europe_22_2020_alignment.sorted.bam.vcf > genbank_sequences_Europe_22_2020_alignment.sorted.bam.DP4


###############################################################
### Final Primer Intersections: Adding USA genbank datasets ###
###############################################################

echo "Final Primer Intersections: Adding USA genbank datasets"
echo ""

cp genbank_sequences_USA_alignment.sorted.bam.vcf ./genbank_USA_March_25_2020.vcf
cp genbank_sequences_Asia_22_2020_alignment.sorted.bam.vcf ./genbank_Asia_April_22_2020.vcf
cp genbank_sequences_Europe_22_2020_alignment.sorted.bam.vcf ./genbank_Europe_April_22_2020.vcf
cp genbank_sequences_North_America_22_2020_alignment.sorted.bam.vcf ./genbank_North_America_April_22_2020.vcf

bgzip genbank_USA_March_25_2020.vcf
tabix -p vcf genbank_USA_March_25_2020.vcf.gz
bgzip genbank_Asia_April_22_2020.vcf
tabix -p vcf genbank_Asia_April_22_2020.vcf.gz
bgzip genbank_Europe_April_22_2020.vcf
tabix -p vcf genbank_Europe_April_22_2020.vcf.gz
bgzip genbank_North_America_April_22_2020.vcf
tabix -p vcf genbank_North_America_April_22_2020.vcf.gz

vcf-merge genbank_USA_March_25_2020.vcf.gz genbank_Asia_April_22_2020.vcf.gz genbank_Europe_April_22_2020.vcf.gz genbank_North_America_April_22_2020.vcf.gz > final_merge.vcf

bedtools intersect -a CDC_primers.bed -b final_merge.vcf > CDC_primers.intersection.genbank
bedtools intersect -a HK_Pasteur_Korea.bed -b final_merge.vcf > HK_Pasteur_Korea.intersection.genbank
bedtools intersect -a primer-blast-50-170-bp.bed -b final_merge.vcf > primer-blast-50-170-bp.intersection.genbank

echo "All Done"
#
} | tee logfile
#

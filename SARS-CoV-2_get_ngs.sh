#!/bin/bash

{

SRA_list=${1}
Reference=${2}
Threads=${3}
path_to_perl5_lib=${4}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads] [path_to_perl5_lib]"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  echo "[path_to_perl5_lib]: Path to PERL5LIB, in VCFtools folder. If vcftools is installed in /home/user/, will be: /home/user/vcftools/src/perl "
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads] [path_to_perl5_lib]"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  echo "[path_to_perl5_lib]: Path to PERL5LIB, in VCFtools folder. If vcftools is installed in /home/user/, will be: /home/user/vcftools/src/perl "
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads] [path_to_perl5_lib]"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  echo "[path_to_perl5_lib]: Path to PERL5LIB, in VCFtools folder. If vcftools is installed in /home/user/, will be: /home/user/vcftools/src/perl "
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads] [path_to_perl5_lib]"
  echo ""
  echo "This program will call variants using SAMtools/bcftools in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  echo "[path_to_perl5_lib]: Path to PERL5LIB, in VCFtools folder. If vcftools is installed in /home/user/, will be: /home/user/vcftools/src/perl "
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads] [path_to_perl5_lib]"; exit 1; }

if [ $# -ne 4c ]; then
  echo 1>&2 "Usage: ./`basename $0` [SRA_list] [Reference] [Threads] [path_to_perl5_lib]"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

echo "Getting Illumina SRA sequencing data"
while read line; do
  fastq-dump -Z --gzip $line > $line.fq.gz
done <${1}
echo" Done"

echo "Trimming reads with fastp"
SRA= ls -1 *.fq.gz
for SRA in *.fq.gz; do fastp -w ${3} -i ${SRA} -o ${SRA}.fastp
done
echo "Done"

echo "Mapping reads againts SARS-CoV-2 reference genome with minimap2"
fastp= ls -1 *.fq.gz.fastp
for fastp in *.fq.gz.fastp; do minimap2 -ax sr ${2} ${fastp} > ${fastp}.sam -t 20
done
echo "Done"

echo "Sorting SAM files, using n threads"
sam= ls -1 *.sam
for sam in *.sam; do samtools sort ${sam} > ${sam}.sorted.bam -@ ${3}
done
echo "Done"

echo "Cleaning intermediate files"
rm *.fastp
rm *.sam
echo "Done"

echo "Renaming files in bash"
for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fq.gz//g')";  done
for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastp.sam//g')";  done
echo "Done"

echo "Indexing bam files"
bam= ls -1 *.bam
for bam in *.bam; do samtools index ${bam} -@ 20
done
echo "Done"

echo "Calling and filtering variants by using bcftools"

bam= ls -1 *.bam
for bam in *.bam; do bcftools mpileup --min-ireads 3 -B -C 50 -d 250 --fasta-ref ${2} --threads ${3} -Ou ${bam}| bcftools call -mv -Ov -o ${bam}.vcf
done

bcf= ls -1 *.sorted.bam.vcf
for bcf in *.sorted.bam.vcf; do bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 2' ${bcf} > ${bcf}.filtered
done
echo "Done"

echo "BGZIP and Tabix founder variants"

founder= ls -1 *.sorted.bam.vcf.filtered
for founder in *.sorted.bam.vcf.filtered; do bgzip ${founder}
done

founder= ls -1 *.sorted.bam.vcf.filtered.gz
for founder in *.sorted.bam.vcf.filtered.gz; do tabix -p vcf ${founder}
done
echo "Done"

echo "Merging founder variants across samples"
export PERL5LIB=${4}
vcf-merge --remove-duplicates --trim-ALTs $(ls -1 *.sorted.bam.vcf.filtered.gz | perl -pe 's/\n/ /g') > founder.vcf
echo "Done"

echo "Summarize genotypes in founder variants" 
vcffixup founder.vcf > founder.fixup.vcf
echo "Done"

echo "Plotting Variants"

R
# install.packages("vcfR")
library(vcfR)
my_vcf <- read.vcfR("founder.fixup.vcf", verbose = FALSE)
chrom <- create.chromR(name="SARS-CoV-2 founder variants", vcf=my_vcf)
chrom <- proc.chromR(chrom, verbose=TRUE)
pdf('vcfR_plot.pdf')
plot(chrom)
dev.off()
pdf('chromoqc_plot.pdf')
chromoqc(chrom, xlim=c(1, 29903))
dev.off()
q()
n
####


echo ""
echo "All done."
echo ""

###############################################################
#
} | tee logfile
#

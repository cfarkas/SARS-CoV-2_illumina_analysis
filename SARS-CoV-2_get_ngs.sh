#!/bin/bash

{

SRA_list=${1}
Reference=${2}
Threads=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

echo "Downloading SRA files from the given list of accessions"
prefetch --max-size 800G -O ./ --option-file ${1}
echo "SRA files were downloaded in current directory"
echo ""
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
for a in *.fastq.gz; do fastp -w ${3} -i ${a} -o ${a}.fastp
done

###########################################################################################
# Aligning illumina datasets againts reference with minimap, using 20 threads (-t option) #
###########################################################################################

echo "Aligning illumina datasets againts reference with minimap, using n threads."
echo ""

b= ls -1 *.fastq.gz.fastp
for b in *.fastq.gz.fastp; do minimap2 -ax sr ${2} ${b} > ${b}.sam -t ${3}
samtools sort ${b}.sam > ${b}.sam.sorted.bam -@ ${3}
rm ${b}.sam
rm ${b}
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

##################################
# Cleaning Up intermediate files #
##################################

echo "Cleaning up intermediate files"
rm *.sam 

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
    --referenceFasta ${2} \
    --runDir ${bam}.founder \
    --exome
# execution on a single local machine with 20 parallel jobs
${bam}.founder/runWorkflow.py -m local -j ${3}
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
echo ""
echo "Filtered founder variants are located in working directory"
echo ""

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

###############################################################
#
} | tee logfile
#

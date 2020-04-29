# SARS-CoV-2_illumina_analysis
Computational analysis to discover founder variants in illumina NGS data from SARS-CoV-2: 

1) This commands will download illumina datasets available in SRA archive corresponding to SARS-CoV-2 untill early April 2020 (please see: https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/) and will obtain founder variants per sample by applying variant calling (by using strelka variant calling). 

2) This founder variants can be also cross-validated by using bcftools (with strict filtering) and snippy variant discovery tool (please see: https://github.com/tseemann/snippy) since the latter pipeline employ bayesian-based variant calling (please see: https://github.com/ekg/freebayes). 

3) We also provide a way to evalute the performance of primer sets currently used for viral testing (see CDC_primers.fasta file) (https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html). 

4) Number of threads were setted to 20 in commands, but it can be increased/decreased. 

# Preeliminars: 

### Installing bowtie2-aligner (for details, see: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
```
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download
unzip download
rm download
cp ./bowtie2-2.3.3.1-linux-x86_64/bowtie2* /usr/local/bin/
```
with sudo privileges
```
sudo apt install bowtie2
```

### Installing minimap2 aligner (for install details, please see: https://github.com/lh3/minimap2)
```
### Installing minimap2
git clone https://github.com/lh3/minimap2
# build minimap2
cd minimap2 && make
# with sudo privileges
sudo cp minimap2 /usr/local/bin/
```

### Installing fastp: An ultra-fast all-in-one FASTQ preprocessor (for details, please see: https://github.com/OpenGene/fastp)
```
git clone https://github.com/OpenGene/fastp.git
# build fastp
cd fastp
make
# with sudo privileges
sudo cp fastp /usr/local/bin/
```

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html.
```
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz.1
cd bedtools2
make
# with sudo privileges
sudo cp ./bin/* /usr/local/bin/
```

### Obtaining and Installing VCFtools
Complete instructions can be found in https://vcftools.github.io/downloads.html. Users with privileges can accomplish with sudo as follows: 

```
### Installing vcftools
git clone https://github.com/vcftools/vcftools.git
cd vcftools/
./autogen.sh
export PERL5LIB=/path/to/your/vcftools-directory/src/perl/ 
./configure
make
make install
```

### Installing vcflib
```
git config --global url.https://github.com/.insteadOf git://github.com/
git clone --recursive git://github.com/vcflib/vcflib.git

#Enter vcflib directory and make
cd vcflib
make   # Needs CMake compiler, with sudo privileges do: sudo apt-get install cmake
cp scripts/* /usr/local/bin/
cp bin/* /usr/local/bin/
```

### Obtaining and installing up-to-date SAMtools with htslib (version >= 1.9)
(Old samtools version can also work). Users need to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplished by downloading the three packages, decompressing it, and doing the following:
```
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
bzip2 -d htslib-1.10.2.tar.bz2
tar -xvf htslib-1.10.2.tar
rm htslib-1.10.2.tar
cd htslib-1.10.2    # and similarly for samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools
sudo cp samtools /usr/local/bin/

# Similarly as htslib, samtools and bcftools can be downloaded as follows:

wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
```

Then in a terminal type
>samtools

to check 1.10 version (using htslib v1.10)

### Obtaining SRA toolkit from ncbi (for downloading reads from SRA archive).
```
### Installing SRA toolkit from ncbi
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz
gunzip sratoolkit.2.9.6-ubuntu64.tar.gz
tar -xvf sratoolkit.2.9.6-ubuntu64.tar
cp sratoolkit.2.9.6-ubuntu64/bin/fastq-dump /usr/local/bin/
cp sratoolkit.2.9.6-ubuntu64/bin/prefetch /usr/local/bin/
```

### Obtaining tabix (to bgzip and tabix vcf files).
wget https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
bzip2 -d tabix-0.2.6.tar.bz2
tar -xvf tabix-0.2.6.tar
cd tabix-0.2.6/
cp bgzip tabix /usr/local/bin/


### Obtaining snippy 

Please refer to Torsten Seemann Repo: https://github.com/tseemann/snippy
```
conda install -c conda-forge -c bioconda -c defaults snippy
```
To install miniconda3 in linux:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### Installing jaqcuard
For information, please see :https://jacquard.readthedocs.io/en/v0.42/installation.html
```
pip install jacquard
```

### Installing vcfR library. For documentation: https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html

In R, type

```
install.packages("vcfR")
```

# Quick Start:

To reproduce all computational steps from the paper: https://www.biorxiv.org/content/10.1101/2020.04.09.034462v1 download this repository, provide path to PERL5LIB vcftools folder (located in/src/perl/ in vcftools folder) and execute the given bash script as follows (Control file SRR10971381 must be in this directory to work). 
```
git clone https://github.com/cfarkas/SARS-CoV-2_illumina_analysis.git
cd SARS-CoV-2_illumina_analysis
export PERL5LIB=/path/to/your/vcftools-directory/src/perl/   
./SARS-CoV-2_commands.sh 
```
These lines will execute all the analyses to obtain founder variants, using 20 threads. Users can modify this number in the script by using nano or another text processor. 
The expected output from these commands is: 

- All variants in illumina samples: merged.vcf
- Founder variants in North America: genbank_sequences_North_America_22_2020_alignment.sorted.bam.vcf
- Founder variants in Europe: genbank_sequences_Europe_22_2020_alignment.sorted.bam.vcf
- Founder variants in Asia: genbank_sequences_Asia_22_2020_alignment.sorted.bam.vcf 

To inspect founder variants and plot it against SARS-CoV-2 reference genome, do the following (Using the genbank_sequences_North_America_22_2020_alignment.sorted.bam.vcf file of North America as example):

```
R 
library(vcfR)
my_vcf <- read.vcfR("genbank_sequences_North_America_22_2020_alignment.sorted.bam.vcf", verbose = FALSE)
chrom <- create.chromR(name="SARS-CoV-2 founder variants", vcf=my_vcf)
chrom <- proc.chromR(chrom, verbose=TRUE)
pdf('vcfR_plot.pdf')
plot(chrom)
dev.off()                                                                                                                               pdf('chromoqc_plot.pdf')
chromoqc(chrom, xlim=c(1, 29903))
dev.off()
```
To obtain a phylogenetic tree using variants collected in VCF file, do the following:
```
### vcf2phylip: Convert SNPs in VCF format to PHYLIP, NEXUS, binary NEXUS, or FASTA alignments for phylogenetic analysis

git clone https://github.com/edgardomortiz/vcf2phylip
cp ./vcf2phylip/vcf2phylip.py ./
# MIN_SAMPLES_LOCUS=1
python vcf2phylip.py -i merged.vcf -n -m 1 --fasta
# MIN_SAMPLES_LOCUS=4
python vcf2phylip.py -i merged.vcf -n -m 4 --fasta
```
For more information, please see: https://github.com/edgardomortiz/vcf2phylip

To compare the obtained founder variants using bcftools, do the following for each bam file:
```
# Calling and filtering variants by using bcftools"
echo ""
bam= ls -1 *.bam
for bam in *.bam; do bcftools mpileup --min-ireads 3 -B -C 50 -d 250 --fasta-ref ${2} --threads ${3} -Ou ${bam}| bcftools call -mv -Ov -o ${bam}.vcf
done
# Filtering variants
bcf= ls -1 *.sorted.bam.vcf
for bcf in *.sorted.bam.vcf; do bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${bcf} > ${bcf}.filtered
```

To compare the obtained founder variants using Snippy, do the following (i.e. Using 20 CPUs):

```
# China

snippy --cpus 20 --outdir ./SRR11059940 --ref SARS-CoV-2.gb --se SRR11059940.fastq.gz
snippy --cpus 20 --outdir ./SRR11059944 --ref SARS-CoV-2.gb --se SRR11059944.fastq.gz
snippy --cpus 20 --outdir ./SRR11059945 --ref SARS-CoV-2.gb --se SRR11059945.fastq.gz
snippy --cpus 20 --outdir ./SRR11059946 --ref SARS-CoV-2.gb --se SRR11059946.fastq.gz
snippy --cpus 20 --outdir ./SRR11059947 --ref SARS-CoV-2.gb --se SRR11059947.fastq.gz

# USA
snippy --cpus 20 --outdir ./SRR11241254 --ref SARS-CoV-2.gb --se SRR11241254.fastq.gz
snippy --cpus 20 --outdir ./SRR11241255 --ref SARS-CoV-2.gb --se SRR11241255.fastq.gz
snippy --cpus 20 --outdir ./SRR11247075 --ref SARS-CoV-2.gb --se SRR11247075.fastq.gz
snippy --cpus 20 --outdir ./SRR11247076 --ref SARS-CoV-2.gb --se SRR11247076.fastq.gz
snippy --cpus 20 --outdir ./SRR11247077 --ref SARS-CoV-2.gb --se SRR11247077.fastq.gz
snippy --cpus 20 --outdir ./SRR11247078 --ref SARS-CoV-2.gb --se SRR11247078.fastq.gz
snippy --cpus 20 --outdir ./SRR11278090 --ref SARS-CoV-2.gb --se SRR11278090.fastq.gz
snippy --cpus 20 --outdir ./SRR11278091 --ref SARS-CoV-2.gb --se SRR11278091.fastq.gz
snippy --cpus 20 --outdir ./SRR11278092 --ref SARS-CoV-2.gb --se SRR11278092.fastq.gz
snippy --cpus 20 --outdir ./SRR11278164 --ref SARS-CoV-2.gb --se SRR11278164.fastq.gz
snippy --cpus 20 --outdir ./SRR11278165 --ref SARS-CoV-2.gb --se SRR11278165.fastq.gz
snippy --cpus 20 --outdir ./SRR11278166 --ref SARS-CoV-2.gb --se SRR11278166.fastq.gz
snippy --cpus 20 --outdir ./SRR11278167 --ref SARS-CoV-2.gb --se SRR11278167.fastq.gz
snippy --cpus 20 --outdir ./SRR11278168 --ref SARS-CoV-2.gb --se SRR11278168.fastq.gz

# Australia
snippy --cpus 20 --outdir ./SRR11397714 --ref SARS-CoV-2.gb --se SRR11397714.fastq.gz
snippy --cpus 20 --outdir ./SRR11397715 --ref SARS-CoV-2.gb --se SRR11397715.fastq.gz
snippy --cpus 20 --outdir ./SRR11397716 --ref SARS-CoV-2.gb --se SRR11397716.fastq.gz
snippy --cpus 20 --outdir ./SRR11397717 --ref SARS-CoV-2.gb --se SRR11397717.fastq.gz
snippy --cpus 20 --outdir ./SRR11397718 --ref SARS-CoV-2.gb --se SRR11397718.fastq.gz
snippy --cpus 20 --outdir ./SRR11397719 --ref SARS-CoV-2.gb --se SRR11397719.fastq.gz
snippy --cpus 20 --outdir ./SRR11397720 --ref SARS-CoV-2.gb --se SRR11397720.fastq.gz
snippy --cpus 20 --outdir ./SRR11397721 --ref SARS-CoV-2.gb --se SRR11397721.fastq.gz
snippy --cpus 20 --outdir ./SRR11397728 --ref SARS-CoV-2.gb --se SRR11397728.fastq.gz
snippy --cpus 20 --outdir ./SRR11397729 --ref SARS-CoV-2.gb --se SRR11397729.fastq.gz
snippy --cpus 20 --outdir ./SRR11397730 --ref SARS-CoV-2.gb --se SRR11397730.fastq.gz
```

# Steps for user-provided datasets from SARS-CoV-2

In order to obtain all founder mutations in user-provided SARS-CoV-2 NGS datasets, users need to execute another bash script: SARS-CoV-2_get_ngs.sh, providing: 

- Sequence read archive accessions of each datasets (SRR prefix)
- SARS-CoV-2 fasta reference
- number of threads for calculations 
- path to PERL5LIB from vcftools 

as follows: 
```
git clone https://github.com/cfarkas/SARS-CoV-2_illumina_analysis.git
cd SARS-CoV-2_illumina_analysis
samtools faidx covid19-refseq.fasta
./SARS-CoV-2_get_ngs.sh SRA_list Reference Threads /path/to/your/vcftools-directory/src/perl/
```
For more information about this script, do

```
./SARS-CoV-2_get_ngs.sh -h 
```
(Control file SRR10971381 must be in this directory as well).

To test, we provided a file in this repository called SARS-CoV-2_curated_list_22_04_2020.tabular with updated SARS-CoV-2 next generation sequencing datasets up to April 22, 2020.

Contact: cfarkas@udec.cl, carlosfarkas@gmail.com

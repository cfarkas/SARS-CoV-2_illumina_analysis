# SARS-CoV-2_illumina_analysis
Computational analysis to discover founder variants in illumina NGS data from SARS-CoV-2: 

1) This commands will download illumina datasets available in SRA archive corresponding to SARS-CoV-2 untill early April 2020 (please see: https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/) and will obtain founder variants per sample by applying variant calling (by using bcftools) and strict filtering. 

2) This founder variants can be also cross-validated by using snippy variant discovery (please see: https://github.com/tseemann/snippy) since this pipeline use bayesian-based variant calling (please see: https://github.com/ekg/freebayes). 

3) We also provide a way to evalute the performance of primer sets currently used for viral testing (see CDC_primers.fasta file) (https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html). 

4) Number of threads were setted to 20 in commands, but it can be increased/decreased. A machine with at least 50 GB ram memory is needed for the genbank sequence alignment but number of cores (n=5) can be decreased to n=1. 

# Preeliminars: 

### Obtaining Bowtie2 aligner (for install details, please see: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
```
sudo apt-get update
sudo apt-get install bowtie2
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
git clone https://github.com/vcftools/vcftools.git
./autogen.sh
./configure
make
# with sudo privileges
sudo make install
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
(Old samtools version can also work). Users needs to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplish by downloading the three packages, decompressing it, and doing the following:
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
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz
gunzip sratoolkit.2.9.6-ubuntu64.tar.gz
tar -xvf sratoolkit.2.9.6-ubuntu64.tar
# with sudo privileges
sudo cp sratoolkit.2.9.6-ubuntu64/bin/fastq-dump /usr/local/bin/
```

### Obtaining snippy 

Please refer to Torsten Seemann Repo: https://github.com/tseemann/snippy
```
conda install -c conda-forge -c bioconda -c defaults snippy
```
without sudo privileges:

```
cd $HOME
git clone https://github.com/tseemann/snippy.git
$HOME/snippy/bin/snippy --help
```
### Installing vcfR library. For documentation: https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html

In R, type

```
install.packages("vcfR")
```

# Quick Start:

To reproduce all computational steps from the paper: https://www.biorxiv.org/content/10.1101/2020.04.09.034462v1 dowmload this repository and execute the given bash script as follows:
```
git clone https://github.com/cfarkas/SARS-CoV-2_illumina_analysis.git
cd SARS-CoV-2_illumina_analysis
./SARS-CoV-2_commands.sh 
```
These lines will execute all the analyses. 

# Steps for user-provided datasets from SARS-CoV-2:

In order to obtain all founder mutations in user-provided SARS-CoV-2 NGS datasets, users need to execute another bash script: SARS-CoV-2_get_ngs.sh, providing the Sequence read archive accessions of each datasets (SRR prefix), SARS-CoV-2 fasta reference, number of threads for calculations and path to PERL5LIB from vcftools as follows. 
```
git clone https://github.com/cfarkas/SARS-CoV-2_illumina_analysis.git
cd SARS-CoV-2_illumina_analysis
samtools faidx covid19-refseq.fasta
./SARS-CoV-2_get_ngs.sh SRA_list Reference Threads /path/to/perl5lib 
```
For more information about this script, do

```
./SARS-CoV-2_get_ngs.sh -h 
```

To inspect founder variants and plot it against SARS-CoV-2 reference genome, do the following:
```
R 
library(vcfR)
my_vcf <- read.vcfR("founder.fixup.vcf", verbose = FALSE)
chrom <- create.chromR(name="SARS-CoV-2 founder variants", vcf=my_vcf)
chrom <- proc.chromR(chrom, verbose=TRUE)
pdf('vcfR_plot.pdf')
plot(chrom)
dev.off()                                                                                                                               pdf('chromoqc_plot.pdf')
chromoqc(chrom, xlim=c(1, 29903))
dev.off()
```

Contact: cfarkas@udec.cl, carlosfarkas@gmail.com

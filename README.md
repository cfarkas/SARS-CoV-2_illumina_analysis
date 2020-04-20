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

### Obtaining and installing up-to-date SAMtools, bcftools and htslib (version 1.9)
Old samtools version will not work. Users needs to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplish by downloading the three packages, decompressing it, and doing the following:
```
cd htslib-1.9    # and similarly for bcftools and samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools and bcftools (with sudo privileges)...
sudo cp samtools /usr/local/bin/
```
Then in a terminal type
>samtools<br>bcftools

to check 1.9 versions (using htslib v1.9)

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

# Quick Start:

```
git clone https://github.com/cfarkas/SARS-CoV-2_illumina_analysis.git
cd SARS-CoV-2_illumina_analysis
./SARS-CoV-2_commands.sh 
```
These lines will execute all the analyses. 

Contact: cfarkas@udec.cl, carlosfarkas@gmail.com

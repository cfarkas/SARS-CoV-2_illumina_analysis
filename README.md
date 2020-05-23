# SARS-CoV-2_illumina_analysis
Computational analysis to discover founder variants in illumina NGS data from SARS-CoV-2: 

1) This commands will download illumina datasets available in SRA archive corresponding to SARS-CoV-2 untill early April 2020 (please see: https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/) and will obtain founder variants per sample by applying variant calling (by using strelka variant calling and bcftools). 

2) This founder variants can be also cross-validated by using snippy variant discovery tool (please see: https://github.com/tseemann/snippy) since the latter pipeline employ bayesian-based variant calling (please see: https://github.com/ekg/freebayes). 

3) We also provide a way to evalute the performance of primer sets currently used for viral testing (see CDC_primers.fasta file) (https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html). 

4) Number of threads were setted to 20 in Quick Start commands, but it can be increased/decreased. 

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
```
wget https://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
bzip2 -d tabix-0.2.6.tar.bz2
tar -xvf tabix-0.2.6.tar
cd tabix-0.2.6/
cp bgzip tabix /usr/local/bin/
```

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

### Installing BEDOPS
For information, please see: https://bedops.readthedocs.io/en/latest/content/installation.html#linux
```
wget https://github.com/bedops/bedops/releases/download/v2.4.39/bedops_linux_x86_64-v2.4.39.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2
cp bin/* /usr/local/bin/
```


### Installing vcfR library. For documentation: https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html

In R, type

```
install.packages("vcfR")
```

### Installing MAFFT multi sequence alignment tool
For more information, please see: https://mafft.cbrc.jp/alignment/software/
```
# Install via conda: 
conda install -c bioconda mafft
# with sudo privileges 
sudo apt install mafft
```

### Installing Fasttree: inferring approximately-maximum-likelihood phylogenetic trees
For information, please see: http://www.microbesonline.org/fasttree/
```
# Install via conda: 
conda install -c bioconda fasttree

# Install from Repository (Parallel version used in this repository)
wget http://www.microbesonline.org/fasttree/FastTree
wget http://www.microbesonline.org/fasttree/FastTree.c
gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm
chmod 755 FastTree*
sudo cp FastTreeMP /usr/local/bin/
```


# Quick Start:

To reproduce all computational steps from the paper: https://www.biorxiv.org/content/10.1101/2020.04.09.034462v1 download this repository, provide path to PERL5LIB vcftools folder (located in/src/perl/ in vcftools folder) and execute the given bash script as follows (Control file SRR10971381 must be in this directory to work). 
```
git clone https://github.com/cfarkas/SARS-CoV-2_illumina_analysis.git
cd SARS-CoV-2_illumina_analysis
export PERL5LIB=/path/to/your/vcftools-directory/src/perl/ 
chmod 755 SARS-CoV-2* covid19-refseq.fasta*
./SARS-CoV-2_commands.sh 
```
These lines will execute all the analyses to obtain founder variants, using 20 threads. Users can modify this number in the script by using nano or another text processor. 
The expected output from these commands is: 

- All variants in illumina samples: merged.vcf
- Founder variants in North America: genbank_sequences_North_America_22_2020_alignment.sorted.bam.vcf
- Founder variants in Europe: genbank_sequences_Europe_22_2020_alignment.sorted.bam.vcf
- Founder variants in Asia: genbank_sequences_Asia_22_2020_alignment.sorted.bam.vcf 

## Subsequent analysis

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

To compare the obtained founder variants using bcftools, do the following for each bam file:

```
bam= ls -1 *.bam
for bam in *.bam; do bcftools mpileup -B -C 50 -d 250 --fasta-ref covid19-refseq.fasta --threads 2 -Ou ${bam}| bcftools call -mv -Ov -o ${bam}.vcf
done
### Filtering variants
bcf= ls -1 *.bam.vcf
for bcf in *.bam.vcf; do bcftools filter -e'%QUAL<10 ||(RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || (DP4[0]+DP4[1])/(DP4[2]+DP4[3]) > 0.3' ${bcf} > ${bcf}.filtered
done
```

To compare the obtained founder variants using Snippy, do the following for each trimmed read dataset(i.e. Using 20 CPUs):

```
fastp= ls -1 *.fastq.gz.fastp
for fastp in *.fastq.gz.fastp; do snippy --cpus 20 --outdir ./${fastp}_snippy --ref SARS-CoV-2.gb --se ${fastp}
done
```

To build a Phylogenetic tree of genbank sequences from Asia, Europe and North America, do the following:
```
# Inside SARS-CoV-2_illumina_analysis folder:

# Renaming and join all GenBank sequences
awk '/^>/{print ">Asia_" ++i; next}{print}' < genbank_sequences_Asia_April_22_2020.fasta > asia.fasta
awk '/^>/{print ">Europe_" ++i; next}{print}' < genbank_sequences_Europe_April_22_2020.fasta > europe.fasta
awk '/^>/{print ">North_America_" ++i; next}{print}' < genbank_sequences_North_America_April_22_2020.fasta > north_america.fasta
cat asia.fasta europe.fasta north_america.fasta > all.fasta

# align all.fasta file in the online mafft server: https://mafft.cbrc.jp/alignment/server/large.html or do:
mafft --thread 16 --reorder all.fasta > all.tree # 64 GB RAM needed, decrease number of threads to use less RAM

# Building tree with all.tree output by using fasttree (parallel version)
FastTreeMP -nt all.tree > Fasttree.tree

# Building tree with all.tree output by using fasttree (non-parallel version)
fasttree -nt all.tree > Fasttree_probe

# Phylogenetic tree can be plotted by using the iTOL server: https://itol.embl.de/upload.cgi or FigTree: https://github.com/rambaut/figtree/
```

# Steps for user-provided datasets from SARS-CoV-2

In order to obtain all founder mutations in user-provided SARS-CoV-2 NGS datasets, users need to execute another bash script: SARS-CoV-2_get_ngs.sh, providing: 

- Sequence read archive accessions of each datasets (SRR prefix list, in tabular format)
- SARS-CoV-2 reference in fasta format
- number of threads for calculations 

as follows: 
```
git clone https://github.com/cfarkas/SARS-CoV-2_illumina_analysis.git
cd SARS-CoV-2_illumina_analysis
samtools faidx covid19-refseq.fasta
chmod 755 SARS-CoV-2* covid19-refseq.fasta*
./SARS-CoV-2_get_ngs.sh SRA_list Reference Threads
```
For more information about this script, do

```
./SARS-CoV-2_get_ngs.sh -h 
```
(Control file SRR10971381 must be in this directory as well).

To test, we provided a file in this repository called SARS-CoV-2_curated_list_22_04_2020.tabular with updated SARS-CoV-2 illumina next generation sequencing datasets up to April 22, 2020. Also, an updated list up to May 11, 2020 is available. 

### VCF files:

- bcftools filtered vcf files containing founder mutations are available in the folder ./bcftools_filtered_vcfs/
- Strelka2 output vcf files including the merged vcf file (merged.vcf) are available in the folder ./Strelka2_vcfs/

### Contact
cfarkas@udec.cl, carlosfarkas@gmail.com

### Notes
- Chinese Sample nCoV5 (SRR11059943) largely diverges from SARS-CoV-2 genome and was excluded from our analysis, as also explained here: Shen Z et al., "Genomic diversity of SARS-CoV-2 in Coronavirus Disease 2019 patients.", Clin Infect Dis, 2020 Mar 4. 
- We excluded SRR11426414, SRR11426415, SRR11426416, SRR11426417, SRR11426418 and SRR11426419 (PRJNA615319) since these datasets correspond to recombinant SARS-CoV-2 sequencing instead of clinical samples. 

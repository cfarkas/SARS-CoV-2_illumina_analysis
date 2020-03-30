# SARS-CoV-19_illumina_analysis
Computational analysis to discover germline mutations in illumina NGS data from SARS-CoV-19: 

1) This commands will download illumina datasets available in SRA archive corresponding to SARS-CoV-19 (please see: https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/) and will obtain germline variants per sample by applying variant calling (by using bcftools) and strict filtering. 

2) This variants can be also challenged by using snippy variant discovery (please see: https://github.com/tseemann/snippy) since this pipeline use bayesian-based variant calling (please see: https://github.com/ekg/freebayes). 

3) We also provide a way to evalute the performance of primer sets currently used for viral testing (see CDC_primers.fasta file) (https://www.cdc.gov/coronavirus/2019-ncov/lab/rt-pcr-panel-primer-probes.html). 

# Preeliminars: 

### Obtaining and Installing BEDTools
Complete instructions can be found in https://bedtools.readthedocs.io/en/latest/content/installation.html. Users with privileges can accomplish with sudo: 

```
sudo apt-get install bedtools
```

### Obtaining and installing up-to-date SAMtools, bcftools and htslib (version 1.9)
Old samtools version will not work. Users needs to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplish by downloading the three packages, decompressing it, and doing the following:
```
cd htslib-1.9    # and similarly for bcftools and samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools and bcftools...
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
sudo cp sratoolkit.2.9.6-ubuntu64/bin/fastq-dump /usr/local/bin/
```

### Obtaining snippy 

Please refer to Torsten Seemann Repo: https://github.com/tseemann/snippy
```
conda install -c conda-forge -c bioconda -c defaults snippy
```

# Quick Start:
```
```

In SARS-CoV-19_illumina_analysis folder (after checking all programs referred in preeliminars), open a terminal do the following: 

```
bash SARS-CoV-19_commands 
```
This sentence will execute all the analyses. 

Contact: cfarkas@udec.cl, carlosfarkas@gmail.com

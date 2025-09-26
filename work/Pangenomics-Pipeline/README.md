# Pangenomic Automated Workflow

## Overview
This project provides an automated pipeline for pangenomic analysis of bacterial strains. The workflow processes genomic data to:
- Cluster homologous genes
- Classify genes into core/accessory/unique categories
- Perform functional annotation
- Generate visualizations

## Installation

### Requirements
- Python 3.6+
- Linux environment (tested on Ubuntu 20.04)
  
### Dependencies
Make sure to install these dependencies before proceeding.
```bash
# Python libraries
pip install pandas numpy scipy matplotlib biopython 
````
### External tools
cd-hit: gene clustering tool

eggnog_mapper: gene annotation tool

The download process for both of these tools can be found on https://github.com/weizhongli/cdhit/releases and https://github.com/eggnogdb/eggnog-mapper respectively.

```bash
# Allow the external tools to be accessed by the code
export PATH=$PATH:/path/to/eggnog-mapper
source ~/.bashrc

export PATH=$PATH:/path/to/cd-hit
source ~/.bashrc
````

## Steps

### Your input
In the codes folder exisits all the scripts necessary to run the workflow. The main script to be executed is workflow1.sh. There is also a folder named AcintoAcinetobacter_baumannii_test. Replace the name of this folder with the name of the species you are working on, for example: Citrobacter_freundii. Inside the folder there is an example DNA coding sequence. Replace the file with you CDS fasta file and name it in this format: Citrobacter_feundii_CDS.fasta
Make sure the CDS.fasta file is a concatenated file of all the strains of the species. If you want to run the pipeline for several species, do the same process. 

### CSV files
In the genomes_ids_fasta folder you will find 3 dummy files each for every species being studied. Replace the name of the files with the name of your species in this format Citrobacter_freundii_genome_ids.csv

### Activate and run workflow1.sh
```bash
chmod +x workflow1.sh
./workflow1.sh

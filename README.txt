Cell types aging project (R code). 

This project tests for cell type enrichment in a orded list of genes. Current datasets are setup to rank genes according to how much they change across the human lifespan. It also tests for gene ontology specific groups that are limited to cell type specific genes (GOGroupByCellTypeTests.R). 

The code is setup for multicore machines using the doParallel R package. 

The main entry R script is mainCode.R. It has a number of parameters and can be ran from the command line (for use on a cluster). Two separate scripts are used for loading Tasic and Zeisel single cell datasets (ZeiselProcessingCode.R and TasicProcessingCode.R)

Most of the data is checked into github with the code. Except the following large files:

Download data at:

Allen data instructions from project base folder:
mkdir "data/Allen data"
mkdir "data/Allen data/case study"
cd "data/Allen data/case study"
wget casestudies.brain-map.org/celltax/data/data_download.zip
unzip data_download.zip

Zeisel data instructions from project base folder::
mkdir "data/Zeisel"
cd "data/Zeisel"
wget http://linnarssonlab.org/blobs/cortex/expression_mRNA_17-Aug-2014.txt

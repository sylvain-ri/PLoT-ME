# ReSPLIT-ME 
RefSeq Segmentation for Pre-classification of Long reads : Identification of Taxonomy with high Memory Efficiency <br>
Sylvain Riondet, NUS/SoC, GIS/Biopolis, Singapore


## Description
#### Pre-processing
- Segmentation of RefSeq into clusters
- Building of taxonomic classifiers' indexes for each cluster

#### Processing of mock communities / metagenomics fastq files
- Pre-classification of long DNA reads (Nanopore/PacBio) to each cluster
- Classification by the classifier with a segment of RefSeq

Only Kraken2 is currently automated, but any classifier able to build its index on .fna files with a provided taxid should work

## Take-aways
- High reduction in memory needs, defined by the number of clusters (Mini Batch K-Means)
- Compatible with existing taxonomic classifiers
- Over-head of the pre-classification (currently ~3-5x, improvements for future releases)

## Usage
#### Pre-processing
`python3 /path/to/Reads_Binning/prod/parse_DB.py <path/RefSeq> <folder/for/clusters> <path/taxonomy> 
 -k 4 -w 10000 -n 10 -o <nameOfFolderToOmit>` <br>
For the full help: `parse_DB.py -h`
#### Pre-classification + classification
`python3 /path/to/Reads_Binning/prod/classify.py <folder/for/clusters> <folder/reports> 
 -i <fastq files to preclassify>` <br>
For the full help: `classify.py -h`


### Requirements
**Tested with:**
- Python == 3.7
- biopython==1.73
- configparser==3.5.3
- ete3==3.1.1
- numpy==1.18.1
- pandas==0.23
- scikit-learn==0.22


### Contact
Author: Sylvain Riondet, PhD student at the National University of Singapore, School of Computing <br>
Email: sylvain@u.nus.edu <br>
Lab: Genome Institute of Singapore <br>
Supervisors: Niranjan Nagarajan & Martin Henz <br>
Date started: 2019-02-21 <br>


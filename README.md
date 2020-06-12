# ReSPLIT-ME 
RefSeq Segmentation for Pre-classification of Long reads : Identification of Taxonomy with high Memory Efficiency <br>
Sylvain Riondet, NUS/SoC, GIS/Biopolis, Singapore

### Pre-processing
- Segmentation of RefSeq into clusters
- Building of taxonomic classifiers' indexes for each cluster

### Processing of mock communities / metagenomics fastq files
- Pre-classification of long DNA reads (Nanopore/PacBio) to each cluster
- Classification by the classifier with a segment of RefSeq

### Take-aways
- High reduction in memory needs 
- Number of clusters controlled by Mini Batch K-Means
- Over-head of the pre-classification (currently ~3-5x)

### Contact
Author: Sylvain Riondet, PhD student at the National University of Singapore, School of Computing <br>
Email: sylvain@u.nus.edu <br>
Lab: Genome Institute of Singapore <br>
Supervisors: Niranjan Nagarajan & Martin Henz <br>
Date started: 2019-02-21 <br>


#!/usr/bin/env bash
# Launch each of the Kraken bins for the designated .fastq/.fasta file.
# $1 is the input file to analyse
# $2 is the output folder
# $3 is the number of cores/threads
# Reports are written into path_reports/<today_date>/$2
# Example ./run_each_kraken_part.sh /home/ubuntu/Disks/HDD500/ONT_Mock_Communities/Mock_10000-uniform-bacteria-l1000-q8.fastq mock_10000

cores=$3
folder_date=$(date +"%Y-%m-%d")
path_reports="/home/ubuntu/Disks/SSD500/Reports"

mkdir -p $path_reports/$folder_date/$2

for cluster in {0..9}
do
echo "Loading Kraken for cluster " $cluster
kraken2 --threads $cores --db /home/ubuntu/Disks/SSD500/Segmentation/Kraken_10_clusters_V1/Kraken2_building/$cluster/ --output $path_reports/$folder_date/$2/$cluster.kraken2_out --report $path_reports/$folder_date/$2/$cluster.kraken2_report $1
done


# For mini kraken2
# kraken2 --threads 2 --db Data/Kraken2_DB/minikraken2_v2_8GB_201904_UPDATE/ /home/ubuntu/Disks/HDD500/ONT_isolates/DRR048282.fastq --output /home/ubuntu/Disks/SSD500/Reports/2019-07-12/isolate_dengue_virus/mini.kraken2_out --report /home/ubuntu/Disks/SSD500/Reports/2019-07-12/isolate_dengue_virus/mini.kraken2_report


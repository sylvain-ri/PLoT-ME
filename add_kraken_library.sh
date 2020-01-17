#!/usr/bin/env bash
# Add the genome segments to kraken libraries before building the hash table

cores=2
for cluster in {0..9}
do
find ~/Disks/SSD500/Segmentation/Kraken_10_clusters_V1/$cluster/*/ -name '*.fa' -print0 | xargs -0 -I{} -n1 -P$cores kraken2-build --add-to-library {} --db /home/ubuntu/Disks/SSD500/Segmentation/Kraken_10_clusters_V1/Kraken2_building/$cluster/
done
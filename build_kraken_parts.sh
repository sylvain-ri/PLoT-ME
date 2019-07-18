cores=2
for cluster in {0..9}
do
kraken2-build --build --threads $cores --db /home/ubuntu/Disks/SSD500/Segmentation/Kraken_10_clusters_V1/Kraken2_building/$cluster/
done
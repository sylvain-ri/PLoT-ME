# Uninstall previous PLoT-Me, build the new one, change the wheel for all linux, install it
python3 -m pip uninstall -y PLoT-ME
python3 -m build
dist_name=$(ls -t dist/*.whl | head -1)
echo "dist_name"
auditwheel repair --plat manylinux2014_x86_64 "dist_name"
wheel_name=$(ls -t wheelhouse/*.whl | head -1)
echo "$wheel_name"
python3 -m pip install "$wheel_name"

path_DB="/mnt/data/PLoT-ME-DB"
rm -r "$path_DB/k4_s10000/"
echo "plot-me.preprocess /mnt/data/NCBI/20190704/refseq /mnt/data/taxonomy $path_DB"

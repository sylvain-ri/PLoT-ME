# Uninstall previous PLoT-Me, build the new one, change the wheel for all linux, install it
python3 -m pip uninstall -y PLoT-ME
python3 -m build
wheel_name=$(find -t dist/*.whl | head -1)
echo "$wheel_name"
auditwheel repair --plat manylinux2014_x86_64 "$wheel_name"
python3 -m pip install "$wheel_name"

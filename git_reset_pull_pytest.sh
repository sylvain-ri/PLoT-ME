# Uninstall previous PLoT-Me, build the new one, change the wheel for all linux, install it
echo " -> GIT RESET"
git reset --hard

echo " -> GIT PULL"
git pull

echo " -> build -e"
python3 -m pip install -e .

echo " -> PYTEST"
python3 -m pytest -v --color=yes --log-level=5

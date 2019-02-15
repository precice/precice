python3 setup.py build_ext --library-dirs=$PRECICE_ROOT/build/last --include-dirs=$PRECICE_ROOT/src
python3 setup.py install --user
python3 setup.py clean --all

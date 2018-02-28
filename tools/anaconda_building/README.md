# Setup and Running

You can use anaconda to Build precice. Helper scripts are provided in this folder.

* install anaconda https://www.anaconda.com/download/#linux
* define ```ANACONDA_ROOT=<path/to/anaconda/installation>``` and ```PRECICE_ROOT=<path/to/precice/folder>```.
* run ```./precice_install.sh```. This file initializes your anaconda environment ```precice```.
* use ```precice_activate.sh``` to activate the environment ```precice``` via ```source precice_activate.sh```.
* run ```python FluidSolver.py precice-config.xml``` and ```python StructureSolver.py precice-config.xml``` each in one shell. Don't forget to activate the environment ```precice_tube```.

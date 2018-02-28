# Setup and Running

You can use anaconda to Build precice. Helper scripts are provided in this folder.

* install anaconda https://www.anaconda.com/download/#linux
* go to ```tools/anaconda_building``` (this folder)
* define ```ANACONDA_ROOT=<path/to/anaconda/installation>``` and ```PRECICE_ROOT=<path/to/precice/folder>``` in ```config.sh```.
* run ```./install.sh```. This file initializes your anaconda environment ```precice```.
* use ```activate.sh``` to activate the environment ```precice``` via ```source activate.sh```.
* run ```./build.sh``` to (clean the potentially existing precice build) and build precice. Output of the cleaning and building process is written to ```scons_clean.log``` and ```scons.log```, respectively.

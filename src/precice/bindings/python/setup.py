import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# check if PRECICE_ROOT is defined
if not os.getenv('PRECICE_ROOT'):
   print("ERROR: PRECICE_ROOT not defined!")
   exit(1)

precice_root = os.getenv('PRECICE_ROOT')


# name of Interfacing API
appname = "PySolverInterface"

# clean previous build
for root, dirs, files in os.walk(".", topdown=False):
   for name in files:
      if ((name.startswith(appname)) and not name.endswith(".pyx")):
         os.remove(os.path.join(root, name))
   for name in dirs:
      if (name == "build"):
         shutil.rmtree(name)

# determine which flags to use with mpic++
if not os.getenv('PRECICE_MPI_IMPLEMENTATION'):
    print('please define PRECICE_MPI_IMPLEMENTATION')
    print('')
    print('use: export PRECICE_MPI_IMPLEMENTATION=<mpi_implementation>')
    print('you can determine your implementation by typing mpic++ -v or mpic++ --showme::version, currently openmpi and mpich are supported.')
    print('')
    print('stopping setup...')
    quit()
elif os.getenv('PRECICE_MPI_IMPLEMENTATION') == 'mpich':
    mpi_compile_args = os.popen("mpic++ -compile-info").read().strip().split(' ')[1::]
    mpi_link_args = os.popen("mpic++ -link-info").read().strip().split(' ')[1::]
elif os.getenv('PRECICE_MPI_IMPLEMENTATION') == 'openmpi':
    mpi_compile_args = os.popen("mpic++ -showme:compile").read().strip().split(' ')
    mpi_link_args = os.popen("mpic++ --showme:link").read().strip().split(' ')
else:
    print('use either mpich or openmpi for PRECICE_MPI_IMPLEMENTATION.')
    print('')
    print('stopping setup...')
    quit()

# need to include libs here, because distutils messes up the order
compile_args = ["-I"+precice_root, "-Wall", "-std=c++11"] + mpi_compile_args
link_args = ["-L"+precice_root+"/build/last/", "-lprecice"] + mpi_link_args

# build precice.so python extension to be added to "PYTHONPATH" later
setup(
   cmdclass = {'build_ext': build_ext},
   ext_modules = [
      Extension(appname,
         sources=[appname+".pyx"],
         libraries=[],
         include_dirs=[precice_root],
         language="c++",
         extra_compile_args=compile_args,
         extra_link_args=link_args
      )
   ]
)

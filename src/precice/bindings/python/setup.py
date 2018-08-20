import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

python_bindings_path = os.path.dirname(os.path.abspath(__file__))
precice_root = os.path.join(python_bindings_path,"../../../../../precice")

# name of Interfacing API
appname = "PySolverInterface"

# determine which flags to use with mpic++
if not os.getenv('PRECICE_MPI_IMPLEMENTATION') or os.getenv('PRECICE_MPI_IMPLEMENTATION') == 'openmpi':  # if PRECICE_MPI_IMPLEMENTATION is not defined try openmpi commands as default (this might break)
    mpi_compile_args = os.popen("mpic++ -showme:compile").read().strip().split(' ')
    mpi_link_args = os.popen("mpic++ --showme:link").read().strip().split(' ')
elif os.getenv('PRECICE_MPI_IMPLEMENTATION') == 'mpich':  # if user provides mpich, then use mpich commands
    mpi_compile_args = os.popen("mpic++ -compile-info").read().strip().split(' ')[1::]
    mpi_link_args = os.popen("mpic++ -link-info").read().strip().split(' ')[1::]
else:  # if PRECICE_MPI_IMPLEMENTATION is not mpich or openmpi quit.
    print('use either mpich or openmpi for PRECICE_MPI_IMPLEMENTATION.')
    print('')
    print('stopping setup...')
    quit()

# need to include libs here, because distutils messes up the order
compile_args = ["-I"+precice_root, "-Wall", "-std=c++11"] + mpi_compile_args
link_args = ["-L"+precice_root+"/build/last/", "-lprecice"] + mpi_link_args

# build precice.so python extension to be added to "PYTHONPATH" later
setup(
    name = appname,
    description = 'Python language bindings for preCICE coupling library',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension(
            appname,
            sources=[os.path.join(python_bindings_path,appname)+".pyx"],
            libraries=[],
            include_dirs=[precice_root],
            language="c++",
            extra_compile_args=compile_args,
            extra_link_args=link_args
        )
    ]
)

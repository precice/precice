import sys
import os
import shutil
import subprocess
from enum import Enum

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

python_bindings_path = os.path.dirname(os.path.abspath(__file__))
precice_root = os.path.join(python_bindings_path,"../../../../../precice")

# name of Interfacing API
appname = "PySolverInterface"


class MpiImplementations(Enum):
    OPENMPI = 1
    MPICH = 2

def check_mpi_implementation():

    FNULL = open(os.devnull, 'w')  # used to supress output of subprocess.call

    if subprocess.call(["mpic++","-showme:compile"], stdout=FNULL, stderr=FNULL) == 0:
        PRECICE_MPI_IMPLEMENTATION = MpiImplementations.OPENMPI
    elif subprocess.call(["mpic++","-compile-info"], stdout=FNULL, stderr=FNULL) == 0:
        PRECICE_MPI_IMPLEMENTATION = MpiImplementations.MPICH
    else:
        raise Exception("unknown/no mpi++")

    return PRECICE_MPI_IMPLEMENTATION


PRECICE_MPI_IMPLEMENTATION = check_mpi_implementation()

# determine which flags to use with mpic++
if PRECICE_MPI_IMPLEMENTATION is MpiImplementations.OPENMPI:
    mpi_compile_args = subprocess.check_output(["mpic++","-showme:compile"]).strip().split(' ')
    mpi_link_args = subprocess.check_output(["mpic++","-showme:link"]).strip().split(' ')
elif PRECICE_MPI_IMPLEMENTATION is MpiImplementations.MPICH:
    mpi_compile_args = subprocess.check_output(["mpic++","-compile-info"]).strip().split(' ')[1::]
    mpi_link_args = subprocess.check_output(["mpic++","-link-info"]).strip().split(' ')[1::]
else:  # if PRECICE_MPI_IMPLEMENTATION is not mpich or openmpi quit.
    raise Exception("unknown/no mpi++. Could not build PySolverInterface.")

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

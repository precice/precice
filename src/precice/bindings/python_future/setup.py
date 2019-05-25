import os
import subprocess
from enum import Enum

from setuptools import setup
from setuptools.command.test import test
import distutils
from distutils.cmd import Command
from Cython.Distutils.extension import Extension
from Cython.Distutils.build_ext import new_build_ext as build_ext
from Cython.Build import cythonize

# name of Interfacing API
APPNAME = "precice_future"
APPVERSION = "1.5.0"  # todo: should be replaced with precice.get_version() as soon as it exists , see https://github.com/precice/precice/issues/261

PYTHON_BINDINGS_PATH = os.path.dirname(os.path.abspath(__file__))

class MpiImplementations(Enum):
    OPENMPI = 1
    MPICH = 2


def check_mpi_implementation(mpi_compiler_wrapper):
    FNULL = open(os.devnull, 'w')  # used to supress output of subprocess.call

    if subprocess.call([mpi_compiler_wrapper, "-showme:compile"], stdout=FNULL, stderr=FNULL) == 0:
        PRECICE_MPI_IMPLEMENTATION = MpiImplementations.OPENMPI
    elif subprocess.call([mpi_compiler_wrapper, "-compile-info"], stdout=FNULL, stderr=FNULL) == 0:
        PRECICE_MPI_IMPLEMENTATION = MpiImplementations.MPICH
    else:
        raise Exception("unknown/no mpi++")

    return PRECICE_MPI_IMPLEMENTATION


def determine_mpi_args(mpi_compiler_wrapper):
    PRECICE_MPI_IMPLEMENTATION = check_mpi_implementation(mpi_compiler_wrapper)
    # determine which flags to use with mpi compiler wrapper
    if PRECICE_MPI_IMPLEMENTATION is MpiImplementations.OPENMPI:
        mpi_compile_args = subprocess.check_output([mpi_compiler_wrapper, "-showme:compile"]).decode().strip().split(
            ' ')
        mpi_link_args = subprocess.check_output([mpi_compiler_wrapper, "-showme:link"]).decode().strip().split(' ')
    elif PRECICE_MPI_IMPLEMENTATION is MpiImplementations.MPICH:
        mpi_compile_args = subprocess.check_output([mpi_compiler_wrapper, "-compile-info"]).decode().strip().split(' ')[
                           1::]
        mpi_link_args = subprocess.check_output([mpi_compiler_wrapper, "-link-info"]).decode().strip().split(' ')[1::]
    else:  # if PRECICE_MPI_IMPLEMENTATION is not mpich or openmpi quit.
        raise Exception("unknown/no mpi found using compiler %s. Could not build PySolverInterface." % mpi_compiler_wrapper)

    return mpi_compile_args, mpi_link_args


def get_extensions(mpi_compiler_wrapper):
    mpi_compile_args, mpi_link_args = determine_mpi_args(mpi_compiler_wrapper)
    
    compile_args = ["-Wall", "-std=c++11"] + mpi_compile_args
    link_args = ["-lprecice"] + mpi_link_args

    return [
        Extension(
                APPNAME,
                sources=[os.path.join(PYTHON_BINDINGS_PATH, APPNAME) + ".pyx"],
                libraries=[],
                language="c++",
                extra_compile_args=compile_args,
                extra_link_args=link_args
            ),
        Extension(
                "test_bindings_module",
                sources=[os.path.join(PYTHON_BINDINGS_PATH, "test", "test_bindings_module" + ".pyx")],
                libraries=[],
                language="c++",
                extra_compile_args=compile_args,
                extra_link_args=link_args
            )
    ]


class my_test(test, object):
    def run(self):
        build_test_package = ['cythonize', '-i', '-E TEST=True', 'precice_future.pyx', 'test/test_bindings_module.pyx']  # before running the tests, we have to build the tests module
        self.announce(
            'Running command: %s' % str(build_test_package),
            level=distutils.log.INFO)
        subprocess.check_call(build_test_package)
        super().run()


# some global definitions for an additional user input command
doc_string = 'specify the mpi compiler wrapper'
mpicompiler_default = "mpic++"

dependencies = ['cython']
dependencies.append('mpi4py')  # only needed, if preCICE was compiled with MPI, see https://github.com/precice/precice/issues/311

setup(
    name=APPNAME,
    version=APPVERSION,
    description='Python language bindings for preCICE coupling library',
    url='https://github.com/precice/precice',
    author='the preCICE developers',
    author_email='info@precice.org',
    license='LGPL-3.0',
    python_requires='>=3',
    install_requires=dependencies,
    ext_modules=cythonize(get_extensions(mpicompiler_default), compile_time_env={"TEST":False}),
    cmdclass={'test': my_test},
    #ensure pxd-files:
    package_data={ 'my_module': ['*.pxd']},
    include_package_data=True,
    zip_safe=False  #needed because setuptools are used
)

import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# check if PRECICE_ROOT is defined
if not os.getenv('PRECICE_ROOT'):
   print "ERROR: PRECICE_ROOT not defined!"
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

# no mpic++ in distutils. So add those flags
mpi_compile_args = os.popen("mpic++ --showme:compile").read().strip().split(' ')
mpi_link_args = os.popen("mpic++ --showme:link").read().strip().split(' ')

# need to include libs here, because distutils messes up the order
compile_args = ["-I"+precice_root, "-Wall", "-std=c++11"] + mpi_compile_args
link_args = ["-L"+precice_root+"/build/last/", "-lprecice", "-lboost_system", "-lboost_filesystem"] + mpi_link_args

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


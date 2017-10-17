# structure/SConstruct
#
# SCons buildfile for "Structure0815" solver
#
######################################

import os;
import sys;
  
preciceRoot = os.getenv ('PRECICE_ROOT')
if ( preciceRoot == None ):
   print('ERROR: Environment variable PRECICE_ROOT not defined!')
   sys.exit(1)
else:
   print('Using environment variable PRECICE_ROOT =', preciceRoot)


##### Declare build variables and default values
#
cppdefines =  []   
libpath = []
libs = [ 'precice' ]
cpppath = [
   preciceRoot + '/src/'
   ]   
ccflags = []
cxx = 'mpicxx'      # For systems offering mpicxx compiler

libpath.append('/usr/lib/')   
libs.append('python2.7')

##### Read command line arguments
#
buildmode = ARGUMENTS.get('build', 'debug')
if buildmode == 'debug':
   cppdefines.append('STRUCTURE_DEBUG_MODE')
   cppdefines.append('Debug')
   cppdefines.append('Asserts')
   ccflags.append('-g3')
   ccflags.append('-O0')
elif buildmode == 'release':
   ccflags.append('-O3')
else:
   print("ERROR: flag 'buildmode' must be set to either 'debug' or 'release'!")
   sys.exit(1)

#libpath.append (preciceRoot + '/build/' + buildmode + '-nopython/')
libpath.append (preciceRoot + '/build/' + buildmode + '/')


##### Setup build environment and issue builds
#
env = Environment ( 
   CPPDEFINES = cppdefines,  # defines for preprocessor (#define xyz)
   LIBPATH    = libpath,     # path to libraries used
   LIBS       = libs,        # libraries used (without prefix "lib" and suffix ".a"/".so"/...)
   CPPPATH    = cpppath,     # pathes where the preprocessor should look for files
   CCFLAGS    = ccflags,     # flags for the c/c++ compilers
   CXX = cxx,                 # the c++ compiler that should be used
   ENV = os.environ
   )
   
env.Program ( buildmode + '/structure0815', Glob('*.cpp') )

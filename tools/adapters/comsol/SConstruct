# comsol/SConstruct
#
# SCons buildfile for "Comsol" executables
#
##########################################

import os;
import sys;
  
env = Environment()
conf = Configure(env) # For checking libraries

##### Determine precice root path from environment variable
#
preciceRoot = os.getenv ('PRECICE_ROOT')
if ( preciceRoot == None ):
   print 'ERROR: Environment variable PRECICE_ROOT not defined!'
   sys.exit(1)
else:
   print 'Using environment variable PRECICE_ROOT =', preciceRoot

##### Determine whether Python scripting extensions should be enabled
#
usePython = ARGUMENTS.get('python', 'on')
pythonLibPaths = []
pythonLibs = []
pythonCppPaths = []
pythonCppDefines = []
if usePython == 'on':
   # Determine Python library path
   pythonLibPath = os.getenv('PRECICE_PYTHON_LIB_PATH')
   if ((pythonLibPath == None) or (pythonLibPath == "")):
      pythonLibPath = '/usr/lib/'
      print 'Env. var. PRECICE_PYTHON_LIB_PATH not defined or empty, assuming python lib to be at "' + pythonLibPath + '"'  
   else:
      print 'Using env. var. PRECICE_PYTHON_LIB_PATH =', pythonLibPath
   pythonLibPaths.append(pythonLibPath)
   
   # Determine Python library name
   pythonLib = os.getenv('PRECICE_PYTHON_LIB')
   if ((pythonLib == None) or (pythonLib == "")):
      pythonLib = 'python2.7'
      print 'Env. var. PRECICE_PYTHON_LIB not defined or empty, assuming python lib to be "' + pythonLib + '"'  
   else:
      print 'Using env. var. PRECICE_PYTHON_LIB =', pythonLib
   pythonLibs.append(pythonLib)
   
   # Determine Python include path
   pythonIncPath = os.getenv('PRECICE_PYTHON_INC_PATH')
   if ((pythonIncPath == None) or (pythonIncPath == "")):
      pythonIncPath = '/usr/include/python2.7/'
      print 'Env. var. PRECICE_PYTHON_INC_PATH not defined or empty, assuming python includes to be at "' + pythonIncPath + '"'  
   else:
      print 'Using env. var. PRECICE_PYTHON_INC_PATH =', pythonIncPath
   pythonCppPaths.append(pythonIncPath)
   
   # Determine NumPy include path
   numpyIncPath = os.getenv('PRECICE_NUMPY_INC_PATH')
   if ((numpyIncPath == None) or (numpyIncPath == "")):
      numpyIncPath = '/usr/include/python2.7/numpy/'
      print 'Env. var. PRECICE_NUMPY_INC_PATH not defined or empty, assuming numpy includes to be at "' + numpyIncPath + '"'  
   else:
      print 'Using env. var. PRECICE_NUMPY_INC_PATH =', numpyIncPath
   pythonCppPaths.append(numpyIncPath) 
elif usePython == 'off':
   buildpath += "-nopython"
   pythonCppDefines.append('PRECICE_NO_PYTHON')
else:
   print "ERROR: argument 'python' must be 'on' or 'off'!"
   sys.exit(1)    

##### Set comsol root path
#
#comsolRoot = '/import/home/student/comsol33'
comsolRoot = '/import/home/software/comsol34'
print 'Comsol root directory = ', comsolRoot

#externalLibpath = '/work/gatzhamm/external/installed/lib'
externalLibpath = '/home/software/comsol_fsi/lib/'

##### Declare build variables for comsol library
#
cppdefinesComsol = []
libpathComsol = [ 
   comsolRoot + '/lib/glnxa64',
   externalLibpath
   ]
libsComsol = [ 
   'glib-2.0',
   'libfsi_com',
   'fsi_mesh',
   'flscriptext' # comsol library
   ]
cpppathComsol = [
   comsolRoot + '/script/external',       # comsol headers
   #'/work/gatzhamm/external/installed/include/',
   '/home/software/comsol_fsi/include/', 
   '/usr/include/glib-2.0/',
   '/usr/lib/glib-2.0/include/',
   '/usr/lib/x86_64-linux-gnu/glib-2.0/include/',
   '/usr/include/mpich2/'
   ]   
ccflagsComsol = []

##### Declare build variables for adapter comsol-precice
#
cppdefinesAdapter = [
   'PRECICE_USE_MPI',
   'Dim2'
   ]   
cppdefinesAdapter.append(pythonCppDefines)

libpathAdapter = [
   externalLibpath,
   '/usr/lib/',
   pythonLibPaths
   ]   
#libpathAdapter.append(pythonLibPaths)
   
libsAdapter = [
   'precice',
   'libfsi_com',
   'fsi_mesh',
   'mpich',
   #'python2.7',
   'glib-2.0'
   ]
libsAdapter.append(pythonLibs)
   
if conf.CheckLib('mpl'):
	libsAdapter.append('mpl')
   
cpppathAdapter = [
   preciceRoot + '/src',                  # precice headers
   '/usr/include/glib-2.0/',
   '/usr/lib/glib-2.0/include/',
   '/usr/lib/x86_64-linux-gnu/glib-2.0/include/',
   '/usr/include/mpich2/',
   '/home/software/comsol_fsi/include/'
   #'/usr/include/python2.7/',
   #'/usr/include/python2.7/numpy/'
   ]   
cpppathAdapter.append(pythonCppPaths)
   
ccflagsAdapter = []


##### Read command line arguments
#
build = ARGUMENTS.get('build', 'debug')
if build == 'debug':
   cppdefinesAdapter.append('Debug')
   cppdefinesAdapter.append('Asserts')
   ccflagsAdapter.append('-g3')
   ccflagsAdapter.append('-O0')
   cppdefinesComsol.append('Debug')
   ccflagsComsol.append('-g3')
   ccflagsComsol.append('-O0')
elif build == 'release':
   ccflagsComsol.append('-O3')
   ccflagsAdapter.append('-O3')
else:
   print "ERROR: flag 'build' must be set to either 'debug' or 'release'!"
   sys.exit(1)

#libpathAdapter.append (peanoSrc + '/build/' + build + '/dim2/lib/')
#libpathAdapter.append (preciceRoot + '/build/tarch/' + build + '/')
libpathAdapter.append (preciceRoot + '/build/' + build + '/')


##### Setup build environment and issue builds
#
envComsol = Environment ( 
   CPPDEFINES = cppdefinesComsol,  # defines for preprocessor (#define xyz)
   LIBPATH    = libpathComsol,     # path to libraries used
   LIBS       = libsComsol,        # libraries used (without prefix "lib" and suffix ".a"/".so"/...)
   CPPPATH    = cpppathComsol,     # pathes where the preprocessor should look for files
   CCFLAGS    = ccflagsComsol,      # flags for the c/c++ compilers
   CC         = 'gcc'
   )
envComsol.SharedLibrary ( build + '/comsol_simulation',
                          'comsol_simulation.c' )
   
envAdapter = Environment (
   CPPDEFINES = cppdefinesAdapter,  # defines for preprocessor (#define xyz)
   LIBPATH    = libpathAdapter,     # path to libraries used
   LIBS       = libsAdapter,        # libraries used (without prefix "lib" and suffix ".a"/".so"/...)
   CPPPATH    = cpppathAdapter,     # pathes where the preprocessor should look for files
   CCFLAGS    = ccflagsAdapter,     # flags for the c/c++ compilers
   CXX = 'g++'             # the c++ compiler that should be used
   )
envAdapterLink = envAdapter.Clone ( CC='g++' )


envAdapter.Object ( build + '/ComsolPrecice', 
                    'ComsolPrecice.c' )
envAdapterLink.Program ( build + '/ComsolPrecice',
                         build + '/ComsolPrecice' )
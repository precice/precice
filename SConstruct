# preCICE/SConstruct 

# Main buildfile for Linux based systems.

import os
import sys

##### Initialize build variables
#
cpppath = [ '#src' ]
cppdefines = [ 'tarch=tarchp2' ]
ccflags = [ '-fPIC' ]
libpath = []
libs = []
linkerflags = []
cxx = 'g++'

env = Environment()
conf = Configure(env) # For checking libraries

##### Determine build directory
#
buildDir = ARGUMENTS.get('builddir', 'build')
buildpath = buildDir + '/'

##### Determine build mode
#
build = ARGUMENTS.get('build', 'debug')
if build == 'debug':
   cppdefines.append('Debug')
   cppdefines.append('Asserts')
   ccflags.append('-g3')
   ccflags.append('-O0')
   buildpath += "debug"
elif build == 'release':
   ccflags.append('-O3')
   buildpath += "release"    
else:
   print "ERROR: Argument 'build' must be either 'debug' or 'release'!"
   sys.exit(1)

print '\nEnvironment variables:'
print '(have to be defined by the user to configure build)'

##### Determine Boost root path from environment variable
#
boostRootPath = os.getenv('BOOST_ROOT')
if (boostRootPath == None):
   print 'ERROR: BOOST_ROOT must be defined!'
   sys.exit(1)
else:
   print 'BOOST_ROOT =', boostRootPath
   cpppath.append(boostRootPath)
   
##### Determine path to peano root from environment variable
#
tarchSrc = os.getenv('TARCH_SRC')
if ((tarchSrc == None) or (tarchSrc == "")):
   print 'TARCH_SRC = ./src/ (default)'
   tarchSrc = './src/'
else:
   print 'TARCH_SRC =', tarchSrc
   cpppath.append(tarchSrc)

##### Determine, whether MPI should be used
#
useMPI = ARGUMENTS.get('mpi', 'on')
if useMPI == 'on':
   # Determine MPI library path
   mpiLibPath = os.getenv('PRECICE_MPI_LIB_PATH')
   if ((mpiLibPath == None) or (mpiLibPath == "")):
      mpiLibPath = '/usr/lib/'
      print 'PRECICE_MPI_LIB_PATH = ' + mpiLibPath + ' (default)'  
   else:
      print 'PRECICE_MPI_LIB_PATH =', mpiLibPath
   libpath.append(mpiLibPath)
   
   # Determine MPI library name
   mpiLib = os.getenv('PRECICE_MPI_LIB')
   if ((mpiLib == None) or (mpiLib == "")):
      mpiLib = 'mpich'
      print 'PRECICE_MPI_LIB = ' + mpiLib + ' (default)'
   else:
      print 'PRECICE_MPI_LIB =', mpiLib   
   libs.append(mpiLib)
   if conf.CheckLib('rt'):
      libs.append('rt') # To be compatible with tarch::utils::Watch clock_gettime
   if (mpiLib == 'mpich'): # MPICH1/2 library
      if conf.CheckLib('mpl'):
         libs.append('mpl')          
   elif (mpiLib == 'mpi'): # OpenMPI library
      if conf.CheckLib('mpi_cxx'):
         libs.append('mpi_cxx')
   
   # Determine MPI include path
   mpiIncPath = os.getenv('PRECICE_MPI_INC_PATH')
   if ((mpiIncPath == None) or (mpiIncPath == "")):
      mpiIncPath = '/usr/include/mpich2'
      print 'PRECICE_MPI_INC_PATH = ' + mpiIncPath + ' (default)'
   else:
      print 'PRECICE_MPI_INC_PATH =', mpiIncPath   
   cpppath.append(mpiIncPath) 
elif useMPI == 'off':
   cppdefines.append ('PRECICE_NO_MPI')
   libs.append('rt') # To be compatible with tarch::utils::Watch clock_gettime
   buildpath += "-nompi"
else:
   print "ERROR: Argument 'useMPI' must be either 'on' or 'off'!"
   sys.exit(1)
   
##### Determine use of socket communication
#
useSockets = ARGUMENTS.get('sockets', 'on')
if useSockets == 'on':
    #ccflags.append('-pthread')
    linkerflags.append('-pthread')
elif useSockets == 'off':
    cppdefines.append('PRECICE_NO_SOCKETS')
    buildpath += "-nosockets"
else:
    print "ERROR: Attribute 'sockets' must be either 'on' or 'off'!"
    sys.exit(1)
   
##### Determine activation of SAGA Grid library
#
#useSAGA = ARGUMENTS.get('saga', 'off')
#if useSAGA == 'off':
#    cppdefines.append('PRECICE_NO_SAGA')
#elif useSAGA == 'on':
#    libs.append('saga_package_advert')
#    libs.append('xyz')
#    libpath.append('/opt/saga-1.5.4/lib/')

##### Determine whether Boost.Spirit 2.0 is available
#
hasSpirit2 = ARGUMENTS.get('spirit2', 'on')
if hasSpirit2 == 'on':
   pass
elif hasSpirit2 == 'off':
   cppdefines.append ('PRECICE_NO_SPIRIT2')
   buildpath += "-nospirit2"
else:
   print "ERROR: Argument 'spirit2' must be either 'on' or 'off'!"
   sys.exit(1)

##### Determine whether Python scripting extensions should be enabled
#
usePython = ARGUMENTS.get('python', 'on')
if usePython == 'on':
   # Determine Python library path
   pythonLibPath = os.getenv('PRECICE_PYTHON_LIB_PATH')
   if ((pythonLibPath == None) or (pythonLibPath == "")):
      pythonLibPath = '/usr/lib/'
      print 'PRECICE_PYTHON_LIB_PATH = ' + pythonLibPath + ' (default)'  
   else:
      print 'PRECICE_PYTHON_LIB_PATH =', pythonLibPath
   libpath.append(pythonLibPath)
   
   # Determine Python library name
   pythonLib = os.getenv('PRECICE_PYTHON_LIB')
   if ((pythonLib == None) or (pythonLib == "")):
      pythonLib = 'python2.6'
      print 'PRECICE_PYTHON_LIB = ' + pythonLib + ' (default)'  
   else:
      print 'PRECICE_PYTHON_LIB =', pythonLib
   libs.append(pythonLib)
   
   # Determine Python include path
   pythonIncPath = os.getenv('PRECICE_PYTHON_INC_PATH')
   if ((pythonIncPath == None) or (pythonIncPath == "")):
      pythonIncPath = '/usr/include/python2.6/'
      print 'PRECICE_PYTHON_INC_PATH = ' + pythonIncPath + ' (default)'  
   else:
      print 'PRECICE_PYTHON_INC_PATH =', pythonIncPath
   cpppath.append(pythonIncPath)
   
   # Determine NumPy include path
   numpyIncPath = os.getenv('PRECICE_NUMPY_INC_PATH')
   if ((numpyIncPath == None) or (numpyIncPath == "")):
      numpyIncPath = '/usr/include/python2.6/numpy/'
      print 'PRECICE_NUMPY_INC_PATH = ' + numpyIncPath + ' (default)'  
   else:
      print 'PRECICE_NUMPY_INC_PATH =', numpyIncPath
   cpppath.append(numpyIncPath) 
elif usePython == 'off':
   buildpath += "-nopython"
   cppdefines.append('PRECICE_NO_PYTHON')
else:
   print "ERROR: argument 'python' must be 'on' or 'off'!"
   sys.exit(1)    

##### Determine use of profiling information
#
gprof = ARGUMENTS.get('gprof', 'off') # Read command line parameter
if gprof == 'off':
   pass
elif gprof == 'on':
   ccflags.append('-p')
   ccflags.append('-pg')
   linkerflags.append('-p')
   linkerflags.append('-pg')
   buildpath += "-gprof"
else:
   print "ERROR: Attribute 'gprof' must be = 'on' or 'off'!"
   sys.exit(1)
   
##### Determine activation of statistics computation
#
#computeStatistics = ARGUMENTS.get('statistics', 'off')
#if computeStatistics == 'off':
#   pass
#elif computeStatistics == 'on':
#   cppdefines.append('PRECICE_STATISTICS')
#   buildpath += "-stat"
#else:
#   print "ERROR: Attribute 'statistics' must be = 'on' or 'off'!"
#   sys.exit(1)

##### Determine, which compiler should be used
#
cxx = ARGUMENTS.get('compiler', 'g++')
if cxx=='icc':
   libpath.append ('/usr/lib/')
   libs.append ('stdc++')
   if build == 'debug':
      #      ccflags.append('-Weffc++')
      ccflags.append('-align')   
   elif build == 'release':
      ccflags.append('-w')
      #      ccflags.append('-vec-report') # Gibt aus wenn vectorisiert wurde
      ccflags.append('-fast')
      ccflags.append('-align')
      ccflags.append('-ansi-alias')
elif cxx == 'g++':
   pass
else:
   print "ERROR: Argument 'cxx' must be either 'g++' or 'icc'!"
   sys.exit(1)

##### Determine build path
#
#buildDir = ARGUMENTS.get('builddir', 'build')
buildpathTarch = buildDir + '/tarch/' + build
#buildpath = buildDir + '/' + build + buildpathParallel + buildpathPython

#libpath.append ('#' + buildpathTarch)
libpath.append ('#' + buildpath), 'sourcesParallelDelta'

env = conf.Finish() # Used to check libraries

##### Setup construction environment
#
env = Environment ( 
   CPPDEFINES = cppdefines, # defines for preprocessor (#define xyz)
   LIBPATH    = libpath,    # path to libraries used
   LIBS       = libs,       # libraries used (without prefix "lib" and suffix ".a"/".so"/...)
   CPPPATH    = cpppath,    # pathes where the preprocessor should look for files
   CCFLAGS    = ccflags,    # flags for the c/c++ compilers
   LINKFLAGS  = linkerflags,# flags given to the linker
   CXX        = cxx,        # the c++ compiler that should be used
   ENV        = os.environ  # propagates environment variables to scons  
   )
    
(sourcesTarch) = SConscript (
   tarchSrc + '/tarch/SConscript-preCICE',
   variant_dir = buildpathTarch, 
   duplicate   = 0
   )

(sourcesPreCICE, sourcesPreCICEMain) = SConscript (
   'src/SConscript-linux', 
   variant_dir = buildpath, 
   duplicate   = 0
   )

sourcesBoost = []
if useSockets == 'on':
    print "\nCopy boost sources for socket communication to build ..."
    if not os.path.exists(buildpath + "/boost/"):
       Execute(Mkdir(buildpath + "/boost/"))
    for file in Glob(boostRootPath + "/libs/system/src/*"):
       Execute(Copy(buildpath + "/boost/", file))   
    for file in Glob(boostRootPath + "/libs/thread/src/pthread/*"):
       Execute(Copy(buildpath + "/boost/", file))  
    sourcesBoost = Glob(buildpath + '/boost/*.cpp')
    

lib = env.StaticLibrary (
   target = buildpath + '/libprecice',
   source = [ 
      sourcesTarch, 
      sourcesPreCICE,
      sourcesBoost
      ]
   )
   
bin = env.Program ( 
   target = buildpath + '/binprecice',
   source = [ 
      sourcesTarch, 
      sourcesPreCICEMain,
      sourcesBoost
      ]
   )

Default(lib, bin)
givenBuildTargets = map(str, BUILD_TARGETS)
#print "Build targets before conversion:", buildtargets
for i in range(len(givenBuildTargets)):
   if givenBuildTargets[i] == "lib":
      BUILD_TARGETS[i] = lib[0]
   elif givenBuildTargets[i] == "bin":
      BUILD_TARGETS[i] = bin[0]
      

buildTargets = ""
for target in map(str, BUILD_TARGETS):
   if buildTargets != "":
      buildTargets += ", "
   buildTargets += target


##### Print build summary
#
print
print "Build summary:"
print "  Targets:   " + buildTargets
print "  Options:   builddir   = " + str(buildDir)
print "             compiler   = " + cxx
print "             build      = " + str(build)
print "             mpi        = " + useMPI
print "             sockets    = " + useSockets
#print "             saga       = " + useSAGA
print "             python     = " + usePython 
print "             spirit2    = " + hasSpirit2
print "             gprof      = " + gprof
#print "             statistics = " + computeStatistics
print

print "  Buildpath: " + buildpath
print "  Buildpath tarch: " + buildpathTarch
print
   
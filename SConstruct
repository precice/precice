# preCICE/SConstruct 

# Main buildfile for Linux based systems.

import os
#import sys

##################################################################### FUNCTIONS

# Checks for a library and appends it to env only if not already appended
#
def uniqueCheckLib(conf, lib):
	if conf.CheckLib(lib, autoadd=0):
		conf.env.AppendUnique(LIBS = [lib])
		return True
	else:
		return False
		
def errorMissingLib(lib, usage):
   print "ERROR: Library '" + lib + "' (needed for " + usage + ") not found!"
   Exit(1) 
   
def errorMissingHeader(header, usage):
   print "ERROR: Header '" + header + "' (needed for " + usage + ") not found or does not compile!"
   Exit(1) 

########################################################################## MAIN

##### Initialize build variables
#
#cpppath = []
#cppdefines = []
#ccflags = []
#libpath = []
#libs = []
#linkerflags = []
cxx = 'g++'

env = Environment()   # For configuring build variables
conf = Configure(env) # For checking libraries, headers, ...

env.Append(CPPPATH = ['#src'])
env.Append(CPPDEFINES = ['tarch=tarchp2'])
env.Append(CCFLAGS = ['-fPIC'])

buildDir = ARGUMENTS.get('builddir', 'build')
buildpath = buildDir + '/'


#
#
############### Check build options
#
#
print
print 'Checking build options ...'

build = ARGUMENTS.get('build', 'debug')    
if build != 'debug' and build != 'release': 
   print "ERROR: Option 'build' must be either 'debug' or 'release'! (default: debug)"
   Exit(1)

useMPI = ARGUMENTS.get('mpi', 'on')
if useMPI != 'on' and useMPI != 'off':  
   print "ERROR: Option 'mpi' must be either 'on' or 'off'! (default: on)"
   Exit(1)
   
useSockets = ARGUMENTS.get('sockets', 'on')
if useSockets != 'on' and useSockets != 'off': 
   print "ERROR: Option 'sockets' must be either 'on' or 'off'! (default: on)"
   Exit(1)
   
hasSpirit2 = ARGUMENTS.get('spirit2', 'on')
if hasSpirit2 != 'on' and hasSpirit2 != 'off':
   print "ERROR: Option 'spirit2' must be either 'on' or 'off'! (default: on)"
   Exit(1)
   
usePython = ARGUMENTS.get('python', 'on')
if hasSpirit2 != 'on' and hasSpirit2 != 'off':
   print "ERROR: Option 'python' must be 'on' or 'off'! (default: on)"
   Exit(1)   
   
gprof = ARGUMENTS.get('gprof', 'off') # Read command line parameter
if gprof != 'on' and gprof != 'off':
   print "ERROR: Option 'gprof' must be = 'on' or 'off'! (default: off)"
   Exit(1) 
   
cxx = ARGUMENTS.get('compiler', 'g++')
if cxx != 'g++' and cxx != 'icc':
   print "ERROR: Option 'cxx' must be either 'g++' or 'icc'! (default: g++)"
   Exit(1)
   
print '... done'


#
#
############### Fetch environment variables
#
#

print
print 'Environment variables used for build ...'
print '(have to be defined by the user to configure build)'

boostRootPath = os.getenv('BOOST_ROOT')
if (boostRootPath == None):
   print 'ERROR: BOOST_ROOT must be defined!'
   Exit(1)
else:
   print 'BOOST_ROOT =', boostRootPath
   
tarchSrc = os.getenv('TARCH_SRC')
if ((tarchSrc == None) or (tarchSrc == "")):
   print 'TARCH_SRC = ./src/ (default)'
   tarchSrc = './src/'
else:
   print 'TARCH_SRC =', tarchSrc

if useMPI == 'on':
   mpiLibPath = os.getenv('PRECICE_MPI_LIB_PATH')
   if ((mpiLibPath == None) or (mpiLibPath == "")):
      mpiLibPath = '/usr/lib/'
      print 'PRECICE_MPI_LIB_PATH = ' + mpiLibPath + ' (default)'  
   else:
      print 'PRECICE_MPI_LIB_PATH =', mpiLibPath
   
   # Determine MPI library name
   mpiLib = os.getenv('PRECICE_MPI_LIB')
   if ((mpiLib == None) or (mpiLib == "")):
      mpiLib = 'mpich'
      print 'PRECICE_MPI_LIB = ' + mpiLib + ' (default)'
   else:
      print 'PRECICE_MPI_LIB =', mpiLib   

   mpiIncPath = os.getenv('PRECICE_MPI_INC_PATH')
   if ((mpiIncPath == None) or (mpiIncPath == "")):
      mpiIncPath = '/usr/include/mpich2'
      print 'PRECICE_MPI_INC_PATH = ' + mpiIncPath + ' (default)'
   else:
      print 'PRECICE_MPI_INC_PATH =', mpiIncPath   
   

if useSockets == 'on':
   socketLibPath = os.getenv('PRECICE_SOCKET_LIB_PATH')
   if ((socketLibPath == None) or (socketLibPath == "")):
      socketLibPath = '/usr/lib/'
      print 'PRECICE_SOCKET_LIB_PATH = ' + socketLibPath + ' (default)'  
   else:
      print 'PRECICE_SOCKET_LIB_PATH =', socketLibPath
      
   socketLib = os.getenv('PRECICE_SOCKET_LIB')
   if ((socketLib == None) or (socketLib == "")):
      socketLib = 'pthread'
      print 'PRECICE_SOCKET_LIB = ' + socketLib + ' (default)'
   else:
      print 'PRECICE_SOCKET_LIB =', socketLib
      
   socketIncPath = os.getenv('PRECICE_SOCKET_INC_PATH')
   if ((socketIncPath == None) or (socketIncPath == "")):
      socketIncPath = '/usr/include/'
      print 'PRECICE_SOCKET_INC_PATH = ' + socketIncPath + ' (default)'
   else:
      print 'PRECICE_SOCKET_INC_PATH =', socketIncPath


#useSAGA = ARGUMENTS.get('saga', 'off')
#if useSAGA == 'off':
#    cppdefines.append('PRECICE_NO_SAGA')
#elif useSAGA == 'on':
#    libs.append('saga_package_advert')
#    libs.append('xyz')
#    libpath.append('/opt/saga-1.5.4/lib/')


if usePython == 'on':
   pythonLibPath = os.getenv('PRECICE_PYTHON_LIB_PATH')
   if ((pythonLibPath == None) or (pythonLibPath == "")):
      pythonLibPath = '/usr/lib/'
      print 'PRECICE_PYTHON_LIB_PATH = ' + pythonLibPath + ' (default)'  
   else:
      print 'PRECICE_PYTHON_LIB_PATH =', pythonLibPath
   
   pythonLib = os.getenv('PRECICE_PYTHON_LIB')
   if ((pythonLib == None) or (pythonLib == "")):
      pythonLib = 'python2.7'
      print 'PRECICE_PYTHON_LIB = ' + pythonLib + ' (default)'  
   else:
      print 'PRECICE_PYTHON_LIB =', pythonLib
   
   pythonIncPath = os.getenv('PRECICE_PYTHON_INC_PATH')
   if ((pythonIncPath == None) or (pythonIncPath == "")):
      pythonIncPath = '/usr/include/python2.7/'
      print 'PRECICE_PYTHON_INC_PATH = ' + pythonIncPath + ' (default)'  
   else:
      print 'PRECICE_PYTHON_INC_PATH =', pythonIncPath
   
   # Determine NumPy include path
   numpyIncPath = os.getenv('PRECICE_NUMPY_INC_PATH')
   if ((numpyIncPath == None) or (numpyIncPath == "")):
      numpyIncPath = '/usr/include/python2.7/numpy/'
      print 'PRECICE_NUMPY_INC_PATH = ' + numpyIncPath + ' (default)'  
   else:
      print 'PRECICE_NUMPY_INC_PATH =', numpyIncPath

print '... done'

##### Determine build path
#
buildpathTarch = buildDir + '/tarch/' + build
#libpath.append ('#' + buildpath), 'sourcesParallelDelta'
env.Append(LIBPATH = [('#' + buildpath), 'sourcesParallelDelta'])



#
#
############### Modify environment according to fetched variables
#
#

print
print 'Configuring build variables ...'

env.Replace(ENV = os.environ)


if cxx=='icc':
   env.AppendUnique(LIBPATH = ['/usr/lib/'])
   env.Append(LIBS = ['stdc++'])
   if build == 'debug':
      env.Append(CCFLAGS = ['-align'])
   elif build == 'release':
      env.Append(CCFLAGS = ['-w', '-fast', '-align', '-ansi-alias'])
elif cxx == 'g++':
   pass
env.Replace(CXX = cxx)


if build == 'debug':
   env.Append(CPPDEFINES = ['Debug', 'Asserts'])
   env.Append(CCFLAGS = ['-g3', '-O0'])
   buildpath += "debug"
elif build == 'release':
   env.Append(CCFLAGS = ['-O3'])
   buildpath += "release"    


env.Append(CPPPATH = [boostRootPath])
env.Append(CPPPATH = [tarchSrc])


if useMPI == 'on':
   env.Append(LIBPATH = [mpiLibPath])
   if not uniqueCheckLib(conf, mpiLib):
      errorMissingLib(mpiLib, 'MPI')
   uniqueCheckLib(conf, 'rt') # To work with tarch::utils::Watch::clock_gettime
   if (mpiLib == 'mpich'): # MPICH1/2/3 library
      uniqueCheckLib(conf, 'mpl')
      uniqueCheckLib(conf, 'pthread')
      #conf.CheckLib('pthread')
   elif (mpiLib == 'mpi'): # OpenMPI library
      uniqueCheckLib(conf, 'mpi_cxx')   
   env.Append(CPPPATH = [mpiIncPath])
   if not conf.CheckHeader('mpi.h'):
      errorMissingHeader('mpi.h', 'MPI')
elif useMPI == 'off':
   env.Append(CPPDEFINES = ['PRECICE_NO_MPI'])
   uniqueCheckLib(conf, 'rt') # To work with tarch::utils::Watch::clock_gettime
   buildpath += "-nompi"


if useSockets == 'on':
   env.Append(LIBPATH = [socketLibPath])
   uniqueCheckLib(conf, socketLib)
   env.Append(CPPPATH = [socketIncPath])
   if socketLib == 'pthread':
      if not conf.CheckHeader('pthread.h'):
         errorMissingHeader('pthread.h', 'Sockets')
    #env.Append(LINKFLAGS = ['-pthread']) # Maybe better???
elif useSockets == 'off':
    env.Append(CPPDEFINES = ['PRECICE_NO_SOCKETS'])
    buildpath += "-nosockets"


if hasSpirit2 == 'on':
   pass
elif hasSpirit2 == 'off':
   env.Append(CPPDEFINES = ['PRECICE_NO_SPIRIT2'])
   buildpath += "-nospirit2"


if usePython == 'on':
   env.Append(LIBPATH = [pythonLibPath])
   if not uniqueCheckLib(conf, pythonLib):
      errorMissingLib(pythonLib, 'Python')
   #env.Append(LIBS = [pythonLib])  
   env.Append(CPPPATH = [pythonIncPath, numpyIncPath]) 
   if not conf.CheckHeader('Python.h'):
      errorMissingHeader('Python.h', 'Python')
   if not conf.CheckHeader(['Python.h', 'arrayobject.h']): 
      errorMissingHeader('arrayobject.h', 'Python NumPy')
elif usePython == 'off':
   buildpath += "-nopython"
   env.Append(CPPDEFINES = ['PRECICE_NO_PYTHON'])


if gprof == 'off':
   pass
elif gprof == 'on':
   env.Append(CCFLAGS = ['-p', '-pg'])
   env.Append(LINKFLAGS = ['-p', '-pg'])
   buildpath += "-gprof"
print '... done'

env = conf.Finish() # Used to check libraries
    
    
    
    
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
    print
    print "Copy boost sources for socket communication to build ..."
    if not os.path.exists(buildpath + "/boost/"):
       Execute(Mkdir(buildpath + "/boost/"))
    for file in Glob(boostRootPath + "/libs/system/src/*"):
       Execute(Copy(buildpath + "/boost/", file))   
    for file in Glob(boostRootPath + "/libs/thread/src/pthread/*"):
       Execute(Copy(buildpath + "/boost/", file))  
    sourcesBoost = Glob(buildpath + '/boost/*.cpp')
    print "... done"
    

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
   
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

env = Environment()   # For configuring build variables
conf = Configure(env) # For checking libraries, headers, ...

env.Append(CPPPATH = ['#src'])
env.Append(CPPDEFINES = ['tarch=tarchp2']) # Was (!) needed for linking to Peano 1
env.Append(CCFLAGS = ['-fPIC'])


#---------------------------------------------------------- Check build options

print
print 'Checking build options ...'
print "(to be given on call of scons as 'OPTION=VALUE ...')"
print "OPTION       VALUE (DESCRIPTION)"

buildDir = ARGUMENTS.get('builddir', 'build')
buildpath = buildDir + '/'
print "builddir   = " + str(buildDir) + " (directory holding build files)"

build = ARGUMENTS.get('build', 'debug')    
if build != 'debug' and build != 'release': 
   print "ERROR: Option 'build' must be either 'debug' or 'release'! (default: debug)"
   Exit(1)
print "build      = " + str(build) + " ('debug' or 'release')"

cxx = ARGUMENTS.get('compiler', 'g++')
if cxx != 'g++' and cxx != 'icc' and cxx != 'mpicxx':
   print "ERROR: Option 'compiler' must be either 'g++' or 'icc' or 'mpicxx'! (default: g++)"
   Exit(1)
print "compiler   = " + cxx + " (Compiler used for building. 'g++', 'icc' or 'mpicxx')"

useMPI = ARGUMENTS.get('mpi', 'on')
if useMPI != 'on' and useMPI != 'off':  
   print "ERROR: Option 'mpi' must be either 'on' or 'off'! (default: on)"
   Exit(1)
print "mpi        = " + useMPI + " (Enables MPI-based communication and running coupling tests.)"
if (useMPI == 'off') and (cxx == 'mpicxx'):
   print "ERROR: Option 'compiler' can be set to 'mpicxx' only when using MPI!"
   Exit(1)
   
useSockets = ARGUMENTS.get('sockets', 'on')
if useSockets != 'on' and useSockets != 'off': 
   print "ERROR: Option 'sockets' must be either 'on' or 'off'! (default: on)"
   Exit(1)
print "sockets    = " + useSockets + " (Enables Socket-based communication.)"

useBoostInstallation = ARGUMENTS.get('boost-inst', 'off')
if useBoostInstallation != 'on' and useBoostInstallation != 'off':  
   print "ERROR: Option 'boost-inst' must be either 'on' or 'off'! (default: off)"
   Exit(1)
print "boost-inst = " + useBoostInstallation + " (Enable if Boost is available compiled and installed.)"

useBoostSpirit2 = ARGUMENTS.get('spirit2', 'on')
if useBoostSpirit2 != 'on' and useBoostSpirit2 != 'off':
   print "ERROR: Option 'spirit2' must be either 'on' or 'off'! (default: on)"
   Exit(1)
print "spirit2    = " + useBoostSpirit2 + " (Used for parsing VRML file geometries and checkpointing.)"
   
usePython = ARGUMENTS.get('python', 'on')
if useBoostSpirit2 != 'on' and useBoostSpirit2 != 'off':
   print "ERROR: Option 'python' must be 'on' or 'off'! (default: on)"
   Exit(1)   
print "python     = " + usePython + " (Used for Python scripted solver actions.)"
   
gprof = ARGUMENTS.get('gprof', 'off') # Read command line parameter
if gprof != 'on' and gprof != 'off':
   print "ERROR: Option 'gprof' must be = 'on' or 'off'! (default: off)"
   Exit(1) 
print "gprof      = " + gprof + " (Used in detailed performance analysis.)" 
   
   
print '... done'


#-------------------------------------------------- Fetch environment variables

print
print 'Environment variables used for this build ...'
print '(have to be defined by the user to configure build)'

tarchSrc = os.getenv('PRECICE_TARCH_SRC')
if ((tarchSrc == None) or (tarchSrc == "")):
   print 'PRECICE_TARCH_SRC = ./src/ (default)'
   tarchSrc = './src/'
else:
   print 'PRECICE_TARCH_SRC =', tarchSrc

boostRootPath = os.getenv('PRECICE_BOOST_ROOT')
if ((boostRootPath == None) or (boostRootPath == "")):
   print 'PRECICE_BOOST_ROOT = /usr/include/ (default)'
   boostRootPath = './src/'
else:
   print 'PRECICE_BOOST_ROOT =', boostRootPath

if useBoostInstallation == 'on':
   if useSockets == 'on':
      boostLibPath = os.getenv('PRECICE_BOOST_LIB_PATH')
      if ((boostLibPath == None) or (boostLibPath == "")):
         boostLibPath = '/usr/lib/'
         print 'PRECICE_BOOST_LIB_PATH = ' + boostLibPath + ' (default)'  
      else:
         print 'PRECICE_BOOST_LIB_PATH =', boostLibPath
         
      boostSystemLib = os.getenv('PRECICE_BOOST_SYSTEM_LIB')
      if ((boostSystemLib == None) or (boostSystemLib == "")):
         boostSystemLib = 'boost_system'
         print 'PRECICE_BOOST_SYSTEM_LIB = ' + boostSystemLib + ' (default)'  
      else:
         print 'PRECICE_BOOST_SYSTEM_LIB =', boostSystemLib
      
      boostThreadLib = os.getenv('PRECICE_BOOST_THREAD_LIB')
      if ((boostThreadLib == None) or (boostThreadLib == "")):
         boostThreadLib = 'boost_thread'
         print 'PRECICE_BOOST_THREAD_LIB = ' + boostThreadLib + ' (default)'  
      else:
         print 'PRECICE_BOOST_THREAD_LIB =', boostThreadLib
      
   #boostIncPath = os.getenv('PRECICE_BOOST_INC_PATH')
   #if ((boostIncPath == None) or (boostIncPath == "")):
   #   boostIncPath = '/usr/include/'
   #   print 'PRECICE_BOOST_INC_PATH = ' + boostIncPath + ' (default)'  
   #else:
   #   print 'PRECICE_BOOST_INC_PATH =', boostIncPath

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



#---------------------------- Modify environment according to fetched variables

print
print 'Configuring build variables ...'

env.Replace(ENV = os.environ)

buildpathTarch = buildDir + '/tarch/' + build
#libpath.append ('#' + buildpath), 'sourcesParallelDelta'
env.Append(LIBPATH = [('#' + buildpath)])


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


env.AppendUnique(CPPPATH = [tarchSrc])


if useBoostInstallation == 'on':
   #env.AppendUnique(CPPPATH = [boostIncPath])
   # The socket implementation is based on Boost libs
   if useSockets=='on':
      print boostLibPath
      env.AppendUnique(LIBPATH = [boostLibPath])
      if not uniqueCheckLib(conf, boostSystemLib):
         errorMissingLib(boostSystemLib, 'Boost')
      if not uniqueCheckLib(conf, boostThreadLib):
         errorMissingLib(boostThreadLib, 'Boost')
env.AppendUnique(CPPPATH = [boostRootPath])
if not conf.CheckCXXHeader('boost/array.hpp'):
   errorMissingHeader('boost/array.hpp', 'Boost')
   
   
if useBoostSpirit2 == 'off':
   env.Append(CPPDEFINES = ['PRECICE_NO_SPIRIT2'])
   buildpath += "-nospirit2"
      

if useMPI == 'on':
   if not cxx == 'mpicxx':
      env.AppendUnique(LIBPATH = [mpiLibPath])
      if not uniqueCheckLib(conf, mpiLib):
         errorMissingLib(mpiLib, 'MPI')
      if (mpiLib == 'mpich'): # MPICH1/2/3 library
         uniqueCheckLib(conf, 'mpl')
         uniqueCheckLib(conf, 'pthread')
         #conf.CheckLib('pthread')
      elif (mpiLib == 'mpi'): # OpenMPI library
         uniqueCheckLib(conf, 'mpi_cxx')   
      env.AppendUnique(CPPPATH = [mpiIncPath])
      if not conf.CheckHeader('mpi.h'):
         errorMissingHeader('mpi.h', 'MPI')
elif useMPI == 'off':
   env.Append(CPPDEFINES = ['PRECICE_NO_MPI'])
   buildpath += "-nompi"
uniqueCheckLib(conf, 'rt') # To work with tarch::utils::Watch::clock_gettime


if useSockets == 'on':
   env.AppendUnique(LIBPATH = [socketLibPath])
   uniqueCheckLib(conf, socketLib)
   env.AppendUnique(CPPPATH = [socketIncPath])
   if socketLib == 'pthread':
      if not conf.CheckHeader('pthread.h'):
         errorMissingHeader('pthread.h', 'Sockets')
    #env.Append(LINKFLAGS = ['-pthread']) # Maybe better???
elif useSockets == 'off':
    env.Append(CPPDEFINES = ['PRECICE_NO_SOCKETS'])
    buildpath += "-nosockets"

if usePython == 'on':
   env.AppendUnique(LIBPATH = [pythonLibPath])
   if not uniqueCheckLib(conf, pythonLib):
      errorMissingLib(pythonLib, 'Python')
   env.AppendUnique(CPPPATH = [pythonIncPath, numpyIncPath]) 
   if not conf.CheckHeader('Python.h'):
      errorMissingHeader('Python.h', 'Python')
   # Check for numpy header needs python header first to compile
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
    
    
    
#--------------------------------------------- Define sources and build targets
    
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
if (useSockets == 'on') and (useBoostInstallation == 'off'):
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
print "Targets:   " + buildTargets
print "Buildpath: " + buildpath
#print "  Buildpath tarch: " + buildpathTarch
print
   
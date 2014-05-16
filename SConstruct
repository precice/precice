# preCICE/SConstruct 

# Main buildfile for Linux based systems.

import os

##################################################################### FUNCTIONS

def uniqueCheckLib(conf, lib):
    """ Checks for a library and appends it to env only if not already appended. """
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

    
def oprint(name, value, description):
    """ Pretty prints an option with value and description. """    
    print "{:10} = {:6} ({})".format(name, value, description)

def vprint(name, value, default=True):
    """ Pretty prints an environment variabe with value and modified or not. """
    mod = "(default)" if default else "(modified)"
    print "{:10} {:25} = {}".format(mod, name, value)


def checkset_var(varname, default):
    """ Checks if environment variable is set, use default otherwise and print the value. """    
    var = os.getenv(varname)
    if not var:
        var = default
        vprint(varname, var)
    else:
        vprint(varname, var, False)
    return var


########################################################################## MAIN

env = Environment()   # For configuring build variables
conf = Configure(env) # For checking libraries, headers, ...

env.Append(CPPPATH = ['#src'])
env.Append(CPPDEFINES = ['tarch=tarchp2']) # Was (!) needed for linking to Peano 1

# Produce position independent code for dynamic linking. makes a difference on the m68k, PowerPC and SPARC.
env.Append(CCFLAGS = ['-fPIC'])


#---------------------------------------------------------- Check build options

print
print 'Checking build options ...'
print "(to be given on call of scons as 'OPTION=VALUE ...')"
oprint("OPTION", "VALUE", "DESCRIPTION")

buildDir = ARGUMENTS.get('builddir', 'build')
buildpath = os.path.join(buildDir, "") # Ensures to have a trailing slash
oprint("builddir", buildDir, "directory holding build files")

build = ARGUMENTS.get('build', 'debug')
if build not in ["debug", "release"]:
    print "ERROR: Option 'build' must be either 'debug' or 'release'! (default: debug)"
    Exit(1)
oprint("build", build, "'debug' or 'release'")

cxx = ARGUMENTS.get('compiler', 'g++')
if cxx != 'g++' and cxx != 'icc' and not cxx.startswith('mpic'):
    print "ERROR: Option 'compiler' must be one of 'g++', 'icc', or starting with 'mpic'! (default: g++)"
    Exit(1)
oprint("compiler", cxx, "Compiler used for building. Can be one of 'g++', 'icc', or starting with 'mpic'.")

useMPI = ARGUMENTS.get('mpi', 'on')
if not useMPI in ["on", "off"]:
    print "ERROR: Option 'mpi' must be either 'on' or 'off'! (default: on)"
    Exit(1)
oprint("mpi", useMPI, "Enables MPI-based communication and running coupling tests.")
if (useMPI == 'off') and cxx.startswith('mpic'):
    print "ERROR: Option 'compiler' must be set to an MPI compiler wrapper only when using MPI!"
    Exit(1)
   
useSockets = ARGUMENTS.get('sockets', 'on')
if not useSockets in ["on", "off"]:
    print "ERROR: Option 'sockets' must be either 'on' or 'off'! (default: on)"
    Exit(1)
oprint("sockets", useSockets, "Enables Socket-based communication.")

useBoostInstallation = ARGUMENTS.get('boost-inst', 'off')
if not useBoostInstallation in ["on", "off"]:
    print "ERROR: Option 'boost-inst' must be either 'on' or 'off'! (default: off)"
    Exit(1)
oprint("boost-inst", useBoostInstallation, "Enable if Boost is available compiled and installed.")

useBoostSpirit2 = ARGUMENTS.get('spirit2', 'on')
if useBoostSpirit2 not in ["on", "off"]:
    print "ERROR: Option 'spirit2' must be either 'on' or 'off'! (default: on)"
    Exit(1)
oprint("spirit2", useBoostSpirit2, "Used for parsing VRML file geometries and checkpointing.")
   
usePython = ARGUMENTS.get('python', 'on')
if usePython not in ["on", "off"]:       
    print "ERROR: Option 'python' must be 'on' or 'off'! (default: on)"
    Exit(1)   
oprint("python", usePython, "Used for Python scripted solver actions.")

gprof = ARGUMENTS.get('gprof', 'off') # Read command line parameter
if gprof not in ["on", "off"]:
    print "ERROR: Option 'gprof' must be = 'on' or 'off'! (default: off)"
    Exit(1) 
oprint("gprof", gprof, "Used in detailed performance analysis.")
   
print '... done'


#-------------------------------------------------- Fetch environment variables

print
print 'Environment variables used for this build ...'
print '(have to be defined by the user to configure build)'

tarchSrc = os.getenv('PRECICE_TARCH_SRC')
if not tarchSrc:
    vprint("PRECICE_TARCH_SRC", "./src/")
    tarchSrc = './src/'
else:
    vprint('PRECICE_TARCH_SRC', tarchSrc, False)

boostRootPath = checkset_var('PRECICE_BOOST_ROOT', "./src")

if useBoostInstallation == 'on':
    if useSockets == 'on':
        boostLibPath = checkset_var('PRECICE_BOOST_LIB_PATH', "/usr/lib/")
        boostSystemLib = checkset_var('PRECICE_BOOST_SYSTEM_LIB', "boost_system")
        boostThreadLib = checkset_var('PRECICE_BOOST_THREAD_LIB', "boost_thread")

      
   #boostIncPath = os.getenv('PRECICE_BOOST_INC_PATH')
   #if ((boostIncPath == None) or (boostIncPath == "")):
   #   boostIncPath = '/usr/include/'
   #   print 'PRECICE_BOOST_INC_PATH = ' + boostIncPath + ' (default)'  
   #else:
   #   print 'PRECICE_BOOST_INC_PATH =', boostIncPath

if useMPI == 'on':
    mpiLibPath = checkset_var('PRECICE_MPI_LIB_PATH', "/usr/lib/")
    
    # Determine MPI library name
    mpiLib = checkset_var('PRECICE_MPI_LIB', "mpich")
    mpiIncPath = checkset_var('PRECICE_MPI_INC_PATH', '/usr/include/mpich2')
   

if useSockets == 'on':
    socketLibPath = checkset_var('PRECICE_SOCKET_LIB_PATH', "/usr/lib")
    socketLib = checkset_var('PRECICE_SOCKET_LIB', "pthread")
    socketIncPath =  checkset_var('PRECICE_SOCKET_INC_PATH', '/usr/include/')


#useSAGA = ARGUMENTS.get('saga', 'off')
#if useSAGA == 'off':
#    cppdefines.append('PRECICE_NO_SAGA')
#elif useSAGA == 'on':
#    libs.append('saga_package_advert')
#    libs.append('xyz')
#    libpath.append('/opt/saga-1.5.4/lib/')


if usePython == 'on':
    pythonLibPath = checkset_var('PRECICE_PYTHON_LIB_PATH', '/usr/lib/')
    pythonLib = checkset_var('PRECICE_PYTHON_LIB', "python2.7")
    pythonIncPath = checkset_var('PRECICE_PYTHON_INC_PATH', '/usr/include/python2.7/')
    numpyIncPath = checkset_var('PRECICE_NUMPY_INC_PATH',  '/usr/include/python2.7/numpy/')

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
    if not cxx.startswith('mpic'):
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
    duplicate = 0
)

(sourcesPreCICE, sourcesPreCICEMain) = SConscript (
    'src/SConscript-linux', 
    variant_dir = buildpath, 
    duplicate = 0
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
   

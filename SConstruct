# preCICE/SConstruct

# Main buildfile for Linux based systems.

import os

##################################################################### FUNCTIONS

def uniqueCheckLib(conf, lib):
    """ Checks for a library and appends it to env if not already appended. """
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

def print_options(vars):
    """ Print all build option and if they have been modified from their default value. """
    for opt in vars.options:
        try:
            is_default = vars.args[opt.key] == opt.default
        except KeyError:
            is_default = True
        vprint(opt.key, env[opt.key], is_default, opt.help)

def vprint(name, value, default=True, description = None):
    """ Pretty prints an environment variabe with value and modified or not. """
    mod = "(default)" if default else "(modified)"
    desc = "   " + description if description else ""
    print "{0:10} {1:25} = {2!s:8}{3}".format(mod, name, value, desc)

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

vars = Variables(None, ARGUMENTS)

vars.Add(PathVariable("builddir", "Directory holding build files.", "build", PathVariable.PathAccept))
vars.Add(EnumVariable('build', 'Build type, either release or debug', "debug", allowed_values=('release', 'debug')))
vars.Add("compiler", "Compiler to use.", "g++")
vars.Add(BoolVariable("omp", "Enables OpenMP-based parallelization of point-to-point communication.", True))
vars.Add(BoolVariable("mpi", "Enables MPI-based communication and running coupling tests.", True))
vars.Add(BoolVariable("sockets", "Enables Socket-based communication.", True))
vars.Add(BoolVariable("boost_inst", "Enable if Boost is available compiled and installed.", False))
vars.Add(BoolVariable("spirit2", "Used for parsing VRML file geometries and checkpointing.", True))
vars.Add(BoolVariable("petsc", "Enable use of the Petsc linear algebra library.", True))
vars.Add(BoolVariable("python", "Used for Python scripted solver actions.", True))
vars.Add(BoolVariable("gprof", "Used in detailed performance analysis.", False))


env = Environment(variables = vars, ENV = os.environ)   # For configuring build variables
# env = Environment(ENV = os.environ)
conf = Configure(env) # For checking libraries, headers, ...


Help(vars.GenerateHelpText(env))
env.Append(CPPPATH = ['#src'])
# env.Append(CPPDEFINES = ['tarch=tarchp2']) # Was (!) needed for linking to Peano 1

# Produce position independent code for dynamic linking. makes a difference on the m68k, PowerPC and SPARC.
env.Append(CCFLAGS = ['-fPIC'])


#---------------------------------------------------------- Check build options

print
print "Build options ..."
print_options(vars)

buildpath = os.path.join(env["builddir"], "") # Ensures to have a trailing slash

if not env["mpi"] and env["compiler"].startswith('mpic'):
    print "ERROR: Option 'compiler' must be set to an MPI compiler wrapper only when using MPI!"
    Exit(1)

print '... done'


#-------------------------------------------------- Fetch environment variables

print
print 'Environment variables used for this build ...'
print '(have to be defined by the user to configure build)'


if env["petsc"]:
    PETSC_DIR = checkset_var("PETSC_DIR", "")
    PETSC_ARCH = checkset_var("PETSC_ARCH", "")

if env["boost_inst"]:
    if env["sockets"]:
        boostLibPath = checkset_var('PRECICE_BOOST_LIB_PATH', "/usr/lib/")
        boostSystemLib = checkset_var('PRECICE_BOOST_SYSTEM_LIB', "boost_system-mt")
        boostThreadLib = checkset_var('PRECICE_BOOST_THREAD_LIB', "boost_thread-mt")
        boostChronoLib = checkset_var('PRECICE_BOOST_CHRONO_LIB', "boost_chrono-mt")
else:
    boostRootPath = checkset_var('PRECICE_BOOST_ROOT', "./src")
      
   #boostIncPath = os.getenv('PRECICE_BOOST_INC_PATH')
   #if ((boostIncPath == None) or (boostIncPath == "")):
   #   boostIncPath = '/usr/include/'
   #   print 'PRECICE_BOOST_INC_PATH = ' + boostIncPath + ' (default)'
   #else:
   #   print 'PRECICE_BOOST_INC_PATH =', boostIncPath

if env["mpi"]:
    mpiLibPath = checkset_var('PRECICE_MPI_LIB_PATH', "/usr/lib/")
    
    # Determine MPI library name
    mpiLib = checkset_var('PRECICE_MPI_LIB', "mpich")
    mpiIncPath = checkset_var('PRECICE_MPI_INC_PATH', '/usr/include/mpich2')
   

if env["sockets"]:
    socketLibPath = checkset_var('PRECICE_SOCKET_LIB_PATH', "/usr/lib")
    socketLib = checkset_var('PRECICE_SOCKET_LIB', "pthread")
    socketIncPath =  checkset_var('PRECICE_SOCKET_INC_PATH', '/usr/include')


#useSAGA = ARGUMENTS.get('saga', 'off')
#if useSAGA == 'off':
#    cppdefines.append('PRECICE_NO_SAGA')
#elif useSAGA == 'on':
#    libs.append('saga_package_advert')
#    libs.append('xyz')
#    libpath.append('/opt/saga-1.5.4/lib/')


if env["python"]:
    pythonLibPath = checkset_var('PRECICE_PYTHON_LIB_PATH', '/usr/lib/')
    pythonLib = checkset_var('PRECICE_PYTHON_LIB', "python2.7")
    pythonIncPath = checkset_var('PRECICE_PYTHON_INC_PATH', '/usr/include/python2.7/')
    numpyIncPath = checkset_var('PRECICE_NUMPY_INC_PATH',  '/usr/include/python2.7/numpy/')

print '... done'


#---------------------------- Modify environment according to fetched variables

print
print 'Configuring build variables ...'

env.Append(LIBPATH = [('#' + buildpath)])
env.Append(CCFLAGS= ['-Wall', '-std=c++11'])


if env["compiler"] == 'icc':
    env.AppendUnique(LIBPATH = ['/usr/lib/'])
    env.Append(LIBS = ['stdc++'])
    if env["build"] == 'debug':
        env.Append(CCFLAGS = ['-align'])
    elif env["build"] == 'release':
        env.Append(CCFLAGS = ['-w', '-fast', '-align', '-ansi-alias'])
elif env["compiler"] == 'g++':
    pass
elif env["compiler"] == "clang++":
    env['ENV']['TERM'] = os.environ['TERM'] # colored compile messages from clang
    env.Append(CCFLAGS= ['-Wsign-compare']) # sign-compare not enabled in Wall with clang.

env.Replace(CXX = env["compiler"])
env.Replace(CC = env["compiler"])

if not conf.CheckCXX():
    Exit(-1)

if env["build"] == 'debug':
    env.Append(CPPDEFINES = ['Debug', 'Asserts'])
    env.Append(CCFLAGS = ['-g3', '-O0'])
    buildpath += "debug"
elif env["build"] == 'release':
    env.Append(CCFLAGS = ['-O3'])
    buildpath += "release"

if env["omp"]:
    if env["compiler"] == 'icc':
        env.Append(CCFLAGS = ['-openmp'])
        env.Append(LINKFLAGS = ['-openmp'])
    elif env["compiler"] == 'g++':
        env.Append(CCFLAGS = ['-fopenmp'])
        env.Append(LINKFLAGS = ['-fopenmp'])
    elif env["compiler"] == "clang++":
        env.Append(CCFLAGS = ['-fopenmp'])
        env.Append(LINKFLAGS = ['-fopenmp'])
    elif env["compiler"].startswith('mpic'):
        env.Append(CCFLAGS = ['-openmp'])
        env.Append(LINKFLAGS = ['-openmp'])
else:
    env.Append(CPPDEFINES = ['PRECICE_NO_OMP'])
    buildpath += "-noomp"

if env["petsc"]:
    env.Append(CPPPATH = [os.path.join( PETSC_DIR, "include"),
                          os.path.join( PETSC_DIR, PETSC_ARCH, "include")])
    env.Append(LIBPATH = [os.path.join( PETSC_DIR, PETSC_ARCH, "lib")])
    if not uniqueCheckLib(conf, "petsc"):
        errorMissingLib("petsc", "Petsc")
else:
    env.Append(CPPDEFINES = ['PRECICE_NO_PETSC'])
    buildpath += "-nopetsc"



if env["boost_inst"]:
    #env.AppendUnique(CPPPATH = [boostIncPath])
    # The socket implementation is based on Boost libs
    if env["sockets"]:
        env.AppendUnique(LIBPATH = [boostLibPath])
    if not uniqueCheckLib(conf, boostSystemLib):
        errorMissingLib(boostSystemLib, 'Boost')
    if not uniqueCheckLib(conf, boostThreadLib):
        errorMissingLib(boostThreadLib, 'Boost')
    if not uniqueCheckLib(conf, boostChronoLib):
        errorMissingLib(boostChronoLib, 'Boost')
else:
    env.AppendUnique(CPPPATH = [boostRootPath])
if not conf.CheckCXXHeader('boost/array.hpp'):
    errorMissingHeader('boost/array.hpp', 'Boost')


if not env["spirit2"]:
    env.Append(CPPDEFINES = ['PRECICE_NO_SPIRIT2'])
    buildpath += "-nospirit2"


if env["mpi"]:
    if not env["compiler"].startswith('mpic'):
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
elif not env["mpi"]:
    env.Append(CPPDEFINES = ['PRECICE_NO_MPI'])
    buildpath += "-nompi"
uniqueCheckLib(conf, 'rt') # To work with tarch::utils::Watch::clock_gettime


if env["sockets"]:
    env.AppendUnique(LIBPATH = [socketLibPath])
    uniqueCheckLib(conf, socketLib)
    env.AppendUnique(CPPPATH = [socketIncPath])
    if socketLib == 'pthread':
        if not conf.CheckHeader('pthread.h'):
            errorMissingHeader('pthread.h', 'Sockets')
    #env.Append(LINKFLAGS = ['-pthread']) # Maybe better???
else:
    env.Append(CPPDEFINES = ['PRECICE_NO_SOCKETS'])
    buildpath += "-nosockets"


if env["python"]:
    env.AppendUnique(LIBPATH = [pythonLibPath])
    if not uniqueCheckLib(conf, pythonLib):
        errorMissingLib(pythonLib, 'Python')
    env.AppendUnique(CPPPATH = [pythonIncPath, numpyIncPath])
    if not conf.CheckHeader('Python.h'):
        errorMissingHeader('Python.h', 'Python')
    # Check for numpy header needs python header first to compile
    if not conf.CheckHeader(['Python.h', 'arrayobject.h']):
        errorMissingHeader('arrayobject.h', 'Python NumPy')
else:
    buildpath += "-nopython"
    env.Append(CPPDEFINES = ['PRECICE_NO_PYTHON'])




if env["gprof"]:
    env.Append(CCFLAGS = ['-p', '-pg'])
    env.Append(LINKFLAGS = ['-p', '-pg'])
    buildpath += "-gprof"

print '... done'


env = conf.Finish() # Used to check libraries

#--------------------------------------------- Define sources and build targets

(sourcesPreCICE, sourcesPreCICEMain) = SConscript (
    'src/SConscript-linux',
    variant_dir = buildpath,
    duplicate = 0
)

sourcesBoost = []
if env["sockets"] and not env["boost_inst"]:
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
    source = [sourcesPreCICE,
              sourcesBoost]
)

bin = env.Program (
    target = buildpath + '/binprecice',
    source = [sourcesPreCICEMain,
              sourcesBoost]
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

print
print "Targets:   " + buildTargets
print "Buildpath: " + buildpath
print

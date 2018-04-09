#!/usr/bin/env python3
"""
Compiles and tests preCICE. Can be used with git bisect run or alike.
Return 0 on success, 1 on failure and 125 on compilation failure which tells git bisect to skip that commit (neither mark it as good or bad)

To ensure a clean test it can delete and recreate the ./tests directory.
"""
import argparse, os, shutil, subprocess, sys

def recreate_testdir(keep_test):
    try:
        if not keep_test:
            print("Removing ./tests")
            shutil.rmtree("./tests")

        print("Creating ./tests")
        os.makedirs("./tests")
    except FileNotFoundError:
        pass
    

def run_test(cmd):
        
    os.chdir("./tests")
    print("Running: ", cmd)
    ret_code = subprocess.call(cmd, shell = True)
    os.chdir("..")
    
    if not ret_code == 0:
        sys.exit(1)


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description = "All additional arguments are passed to the testprecice binary.")

parser.add_argument('-c', '--compile', help="Compile preCICE", dest='compile', action='store_true')
parser.add_argument('-r', '--removebuild', help="Remove build/ and .scon* files before compiling", dest='remove_build', action='store_true')
parser.add_argument('-k', '--keeptest', help="Do not remove test directory for earch test run", dest='keep_test', action='store_true')
parser.add_argument('-t', help="Run tests.", dest='run_tests', action="store_true")
parser.add_argument('-j', help="Number of CPUs to compile on", dest='compile_cpus', default=4)
parser.add_argument('-p', help="Number of MPI procs. Setting to 1 means to not use MPI at all. This does not affect the build process.", type=int, dest='mpi_procs', default=4)
parser.add_argument("-s", "--split", help="Redirect output to a process-unique file testout.N", dest="split_output", action="store_true")
parser.add_argument('--logconfig', "-l", help="Log config file", default = "")
parser.add_argument('--root', help="preCICE Root, defaults to $PRECICE_ROOT", default = os.getenv("PRECICE_ROOT"))
parser.add_argument('--mpirun', help="mpirun executable", default = "mpirun")


if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)
    
args, remainder = parser.parse_known_args()

try:
    os.chdir(args.root)
except TypeError:
    print("The preCICE root directory is not defined. Have you set the $PRECICE_ROOT environment variable?")
    sys.exit(1)
except FileNotFoundError:
    print("$PRECICE_ROOT directory does not exist. Please set the $PRECICE_ROOT environment variable to a valid directory.")
    sys.exit(1)

if args.compile:
    if args.remove_build:
        print("Deleting build directory and scons caches...")
        try:
            shutil.rmtree("./build")
            shutil.rmtree("./.sconf_temp")
            os.remove(".sconsign.dblite")
        except FileNotFoundError:
            pass
    COMPILE_CMD = 'scons mpi=on petsc=on compiler="mpicxx" build=debug -j {cpus}'.format(cpus=args.compile_cpus)
    
    if subprocess.call(COMPILE_CMD, shell = True) != 0:
        sys.exit(125) # Cannot compile, 125 means to skip that revision


mpi_cmd = "{mpirun} -n {procs} {output_filename}".format(mpirun = args.mpirun,
                                                         procs = args.mpi_procs,
                                                         output_filename = "--output-filename 'testout'" if args.split_output else "")

if args.run_tests:
    recreate_testdir(args.keep_test)
    if args.logconfig:
        print("Copy", args.logconfig, "to tests/log.conf")
        shutil.copyfile(args.logconfig, "tests/log.conf")

    run_cmd = "{mpi} ../build/last/testprecice --color_output {args}".format(mpi = mpi_cmd, args = " ".join(remainder))
    run_test(run_cmd)

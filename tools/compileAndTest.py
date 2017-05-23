#!/usr/bin/env python3
"""
Compiles and tests preCICE. Can be used with git bisect run or alike.
Return 0 on success, 1 on failure and 125 on compilation failure which tells git bisect to skip that commit (neither mark it as good or bad)

To ensure a clean test it can delete and recreate the ./tests directory.
"""
import argparse, os, shutil, subprocess, sys

def run_test(cmd, keep_test):
    if not keep_test:
        try:
            print("Removing ./tests")
            shutil.rmtree("./tests")
        except FileNotFoundError:
            pass

    try:
        os.makedirs("./tests")
    except FileExistsError:
        pass
    
    os.chdir("./tests")
    print("Running: ", cmd)
    ret_code = subprocess.call(cmd, shell = True)
    os.chdir("..")
    
    if not ret_code == 0:
        sys.exit(1)


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-c', '--compile', help="Compile preCICE", dest='compile', action='store_true')
parser.add_argument('-r', '--removebuild', help="Remove build/ and .scon* files before compiling", dest='remove_build', action='store_true')
parser.add_argument('-k', '--keeptest', help="Do not remove test directory for earch test run", dest='keep_test', action='store_true')
parser.add_argument('-b', help="Run boost tests.", dest='run_boostttests', action="store_true")
parser.add_argument('-u', help="Run unit tests.", dest='run_unit', action="store_true")
parser.add_argument('-i', help="Run integration tests.", dest='run_integration', action="store_true")
parser.add_argument('-j', help="Number of CPUs to compile on", dest='compile_cpus', default=4)
parser.add_argument('-p', help="Number of MPI procs. Setting to 1 means to not use MPI at all. This does not affect the build process.", type=int, dest='mpi_procs', default=4)
parser.add_argument("-s", "--split", help="Redirect output to a process-unique file testout.N", dest="split_output", action="store_true")
parser.add_argument('--unitconfig', help="Configuration to use for unit tests", dest="unit_test_config", default=".ci-test-config.xml")
parser.add_argument('--integrationconfig', help="Configuration to use for integration tests", dest="integration_test_config", default=".ci-integration-test-config.xml")
parser.add_argument('--logconfig', "-l", help="Log config file", default = "")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)
    
args = parser.parse_args()


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


mpi_cmd = "mpirun -n {procs} {output_filename}".format(procs = args.mpi_procs,
                                                       output_filename = "--output-filename 'testout'" if args.split_output else "")


# Boost Tests
if args.run_boostttests:
    run_cmd = "{mpi} ../build/last/testprecice --color_output".format(mpi = mpi_cmd)
    run_test(run_cmd, args.keep_test)

# Tarch Unit Tests
if args.run_unit:
    run_cmd = "{mpi} ../build/last/binprecice test ../{config} ../src {logconfig}".format(mpi = mpi_cmd,
                                                                                          config = args.unit_test_config,
                                                                                          logconfig = args.logconfig)
    run_test(run_cmd, args.keep_test)

# Tarch Integration Tests
if args.run_integration:
    run_cmd = "{mpi} ../build/last/binprecice test ../{config} ../src {logconfig}".format(mpi = mpi_cmd,
                                                                                          config = args.integration_test_config,
                                                                                          logconfig = args.logconfig)
    run_test(run_cmd, args.keep_test)


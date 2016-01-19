#!/usr/bin/env python
"""
Compiles and tests preCICE. Can be used with git bisect run or alike.
Return 0 on success, 1 on failure and 125 on compilation failure which tells git bisect to skip that commit (neither mark it as good or bad)

To ensure a clean test it deletes and recreates the ./tests directory.
"""
import argparse, os, shutil, subprocess, sys

def run_test(cmd):
    print("Running: ", cmd)

    shutil.rmtree("./tests")
    os.makedirs("./tests")
    os.chdir("./tests")
    proc = subprocess.run(cmd, shell = True)
    os.chdir("..")
    
    if not proc.returncode == 0:
        sys.exit(1)


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-c', '--compile', help="Compile preCICE", dest='compile', action='store_true')
parser.add_argument('-u', help="Run unit tests.", dest='run_unit', action="store_true")
parser.add_argument('-i', help="Run integration tests.", dest='run_integration', action="store_true")
parser.add_argument('-j', help="Number of CPUs to compile on", dest='compile_cpus', default=4)
parser.add_argument('-p', help="Number of MPI procs. Setting to 1 means to not use MPI at all. This does not affect the build process.", type=int, dest='mpi_procs', default=4)
args = parser.parse_args()


if args.compile:
    COMPILE_CMD = 'scons boost_inst=true mpi=on petsc=on compiler="mpicxx" build=debug -j {cpus}'.format(cpus=args.compile_cpus)
    proc = subprocess.run(COMPILE_CMD, shell = True)
    
    if not proc.returncode == 0:
        sys.exit(125) # Cannot compile, 125 means to skip that revision


mpi_cmd = "mpirun -n %s" % args.mpi_procs if args.mpi_procs > 1 else ""

if args.run_unit:
    run_cmd = "{mpi} ../build/last/binprecice test ../{config} ../src".format(mpi = mpi_cmd, config=".ci-test-config.xml")
    run_test(run_cmd)

if args.run_integration:
    run_cmd = "{mpi} ../build/last/binprecice test ../{config} ../src".format(mpi = mpi_cmd, config=".ci-integration-test-config.xml")
    run_test(run_cmd)


#!/usr/bin/env python
"""
Compiles and tests preCICE. Can be used with git bisect run or alike.
Return 0 on success, 1 on failure and 125 on compilation failure which tells git bisect to skip that commit (neither mark it as good or bad)

To ensure a clean test it deletes and recreates the ./tests directory.
"""

import argparse, os, shutil, sys
from subprocess import run

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-n', '--nocompile', help="Disable recompilation", dest='compile', action='store_false')
parser.add_argument('-j', help="Number of CPUs to compile on", dest='compile_cpus', default=4)
args = parser.parse_args()


COMPILE_CMD = 'scons boost_inst=true mpi=on petsc=on compiler="mpicxx" build=debug -j {cpus}'.format(cpus=args.compile_cpus)

shutil.rmtree("./tests")
os.makedirs("./tests")

if args.compile:
    proc = run(COMPILE_CMD, shell = True)
    
    if not proc.returncode == 0:
        sys.exit(125) # Cannot compile, 125 means to skip that revision

os.chdir("./tests")

proc = run("mpirun -n 4  ../build/last/binprecice test ../.ci-test-config.xml ../src", shell = True)

os.chdir("..")

if proc.returncode == 0:
    sys.exit(0)
else:
    sys.exit(1)

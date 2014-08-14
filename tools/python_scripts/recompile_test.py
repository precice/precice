#!env python
""" Clears build tree, rebuilds and runs all the test. Breaks if one command returns a return code != 0.
Could be used to automate search using git bisect. """

import argparse, os, subprocess, sys

parser = argparse.ArgumentParser()

parser.add_argument('-n', '--nocompile', help = "Disable scons -c and recompilation",
                    dest = 'compile', action = 'store_false')

args = parser.parse_args()

cwd = os.getcwd()

try:
    os.mkdir("test_temp")
except OSError:
    pass
    
os.chdir("test_temp")

cmd = [
    "scons --directory .. -c" if args.compile else "",
    "scons --directory .. -j 16 compiler=clang++" if args.compile else "",
    "mpirun.mpich2 -n 4 ../build/debug/binprecice test ../integration-test-config.xml ../src",
    "mpirun.mpich2 -n 4 ../build/debug/binprecice test ../test-config.xml ../src"
]

for c in cmd:
    ret_code = subprocess.call(c, shell=True)
    if ret_code != 0:
        break

os.chdir(cwd)
sys.exit(ret_code)


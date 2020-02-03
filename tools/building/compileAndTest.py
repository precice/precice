#!/usr/bin/env python3
"""
Compiles and tests preCICE. Can be used with git bisect run or alike.
Return 0 on success, 1 on failure and 125 on compilation failure which tells git bisect to skip that commit (neither mark it as good or bad)

To ensure a clean test it can delete and recreate the ./build/TestOutput/ directories.
"""
import argparse, os, shutil, subprocess, sys, multiprocessing


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description = "All additional arguments are passed to the testprecice binary.")

parser.add_argument('-s', "--source-dir", help="Source directory of preCICE", dest='source_dir', default=os.getcwd())
parser.add_argument('-b', "--build-dir", help="Build directory of preCICE", dest='build_dir', default=os.path.join(os.getcwd(), "build"))

parser.add_argument('-c', '--compile', help="Compile preCICE", dest='compile', action='store_true')
parser.add_argument('-r', '--removebuild', help="Remove build/ and .scon* files before compiling", dest='remove_build', action='store_true')
parser.add_argument('-k', '--keeptest', help="Do not remove test directory for each test run", dest='keep_test', action='store_true')
parser.add_argument('-t', help="Run tests.", dest='run_tests', action="store_true")
parser.add_argument('-j', help="Number of CPUs to compile on", dest='compile_cpus', default=multiprocessing.cpu_count())
parser.add_argument('--logconfig', "-l", help="Log config file", default = "")


def wipedir(dir):
    try:
        shutil.rmtree(dir)
        os.makedirs(dir)
    except FileNotFoundError:
        pass

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args, remainder = parser.parse_known_args()


if args.compile:
    if args.remove_build:
        print("Wiping build directory ...")
        wipedir(args.build_dir)
    else:
        os.makedirs(args.build_dir, exist_ok=True)

    CONFIGURE_CMD = 'cmake -DPRECICE_MPICommunication=ON -DPRECICE_PETScMapping=ON -DCMAKE_BUILD_TYPE=Debug -DBUILD_SHARED_LIBS=ON {src}'.format(src=args.source_dir)
    if subprocess.call(CONFIGURE_CMD, shell = True, cwd=args.build_dir) != 0:
        sys.exit(125) # Cannot compile, 125 means to skip that revision

    COMPILE_CMD = 'cmake --build {dir} -- -j {cpus}'.format(dir=args.build_dir, cpus=args.compile_cpus)
    if subprocess.call(COMPILE_CMD, shell = True) != 0:
        sys.exit(125) # Cannot compile, 125 means to skip that revision

if args.run_tests:
    testbase = os.path.join(args.build_dir, "TestOutput")
    testdirs = [os.path.join(testbase, dir) for dir in next(os.walk(testbase))[1]]
    if not args.keep_test:
        print("Wiping test directories ...")
        for dir in testdirs:
            print("Wiping directory: {}".format(dir))
            wipedir(dir)

    if args.logconfig:
        print("Copy ", args.logconfig, " to test directories")
        for dst in testdirs:
            shutil.copyfile(args.logconfig, os.path.join(dir, "log.conf"))

    print("Running Tests:", )
    ret_code = subprocess.call("ctest -VV --output-on-failure", shell = True, cwd=args.build_dir)

    if not ret_code == 0:
        sys.exit(1)

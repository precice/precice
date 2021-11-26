#!/usr/bin/env python3

import os
import argparse
import pathlib


def is_precice_root(dir):
    detect = [
        "CHANGELOG.md", "CMakeLists.txt", "LICENSE", "src", "tests", "cmake"
    ]
    return all(map(lambda c: os.path.exists(os.path.join(dir, c)), detect))


def find_precice_root():
    search_depth = 4
    current = pathlib.Path(os.path.curdir).absolute()
    candidats = [current] + list(current.parents)[:search_depth]
    for dir in candidats:
        if is_precice_root(dir):
            return dir
    raise "Unable to find the root directory of precice"


def patharg(arg):
    assert (arg)
    if " " in arg:
        raise argparse.ArgumentTypeError(
            "The given suite name \"{}\" cannot contain spaces.".format(arg))
    if "." in arg:
        raise argparse.ArgumentTypeError(
            "The given suite name \"{}\" cannot contain the file extensions.".
            format(arg))
    if not arg[0].isupper():
        raise argparse.ArgumentTypeError(
            "The given suite name \"{}\" must use CamelCase.".format(arg))
    return arg


def testname(arg):
    assert (arg)
    if " " in arg:
        raise argparse.ArgumentTypeError(
            "The given test name \"{}\" cannot contain spaces.".format(arg))
    if "." in arg:
        raise argparse.ArgumentTypeError(
            "The given test name \"{}\" cannot contain the file extensions.".
            format(arg))
    if not arg[0].isupper():
        raise argparse.ArgumentTypeError(
            "The given test name \"{}\" must use CamelCase.".format(arg))
    return arg


PRECICE_TEST_BODY = """
{
  PRECICE_TEST(TODO);
  precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);

  std::vector<VertexID> vertexIDs;

  if (context.isNamed(TODO)) {
    auto meshID = interface.getMeshID(TODO);
    auto dataID = interface.getDataID(TODO, meshID);

    std::vector<double> coords;

    interface.setMeshVertices(meshID, TODO, coords.data(), vertexIDs.data());
  } else {
  }
}
"""


def generateTestSource(name, suite, filepath):
    if os.path.exists(filepath):
        raise "The test source at \"{}\" already exists.".format(filepath)

    includes = [
        "<precice/SolverInterface.hpp>", "<vector>", '"testing/Testing.hpp"'
    ]
    suites = ["PreciceTests"] + suite
    space = [""]
    lines = ["#ifndef PRECICE_NO_MPI"]
    lines += space
    lines += ["#include " + inc for inc in includes if inc[0] != "<"]
    lines += space
    lines += ["#include " + inc for inc in includes if inc[0] == "<"]
    lines += space
    lines += ["BOOST_AUTO_TEST_SUITE({})".format(s) for s in suites]
    lines += ["BOOST_AUTO_TEST_CASE({})".format(name), PRECICE_TEST_BODY]
    lines += ["BOOST_AUTO_TEST_SUITE_END() // " + s for s in suites]
    lines += space
    lines += ["#endif // PRECICE_NO_MPI"]

    with open(filepath, "w") as f:
        f.writelines([line + "\n" for line in lines])


def generateTestConfig(name, suite, filepath):
    if os.path.exists(filepath):
        print("The test config at \"{}\" already exists.".format(filepath))
    else:
        open(filepath, 'w').close()


def main():
    parser = argparse.ArgumentParser(
        description="preCICE integration test creation tool.")
    parser.add_argument("suite",
                        metavar="TestSuite",
                        nargs="+",
                        type=patharg,
                        help="test suite(s) to add the test to")
    parser.add_argument("name",
                        metavar="TestName",
                        type=testname,
                        help="name of the test")
    parser.add_argument("-n",
                        "--dry-run",
                        action="store_true",
                        help="print actions only")
    args = parser.parse_args()

    root = find_precice_root()
    location = os.path.join(root, "tests", *args.suite)

    print("Create directory {}".format(location))
    if not args.dry_run:
        os.makedirs(location, exist_ok=True)

    source = args.name + ".cpp"
    config = args.name + ".xml"
    sourcePath = os.path.join(location, source)
    configPath = os.path.join(location, config)

    print("Create test source {}".format(source))
    if not args.dry_run:
        generateTestSource(args.name, args.suite, sourcePath)

    print("Create test config {}".format(config))
    if not args.dry_run:
        generateTestConfig(args.name, args.suite, configPath)


if __name__ == '__main__':
    main()

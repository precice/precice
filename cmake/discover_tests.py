import argparse
import pathlib
import subprocess
import re


def get_tests():
    out = subprocess.run(
        [args.executable, "--list_units"], stdout=subprocess.PIPE, check=True
    ).stdout
    for line in out.strip().splitlines():
        name, ranks = line.decode().split()
        yield name, int(ranks)


def test_dir(root: pathlib.Path, testname: str):
    return root / re.sub(r"[\[\]<>/]", "_", testname)


def get_labels(testname: str):
    suite = testname[: testname.find("/")].removesuffix("Tests")

    labels = [suite] + [
        kw
        for kw in [
            "Ginkgo",
            "Sockets",
            "MPIPorts",
            "MPISinglePorts",
            "MPI",
            "Remeshing",
        ]
        if kw in testname
    ]

    return list(map(str.lower, labels))


def add_command(out, cmd: str, cargs: list[str]):
    print(f"{cmd}({' '.join(cargs)})", file=out)


parser = argparse.ArgumentParser()
parser.add_argument("--executable", type=pathlib.Path)
parser.add_argument("--output", type=argparse.FileType("w"))
parser.add_argument("--run-dir", type=pathlib.Path)
parser.add_argument("--mpi", action="store_true")
parser.add_argument("--mpi-version")
parser.add_argument("--mpi-exec", type=pathlib.Path)
parser.add_argument("--mpi-extra")
parser.add_argument("--mpi-np")
parser.add_argument("--mpi-pre")
parser.add_argument("--mpi-post")
parser.add_argument("--paths")

args = parser.parse_args()

for name, ranks in get_tests():

    print(f"# {name} {ranks}", file=args.output)

    if (
        args.mpi
        and ("MPIPorts" in name or "MPISinglePorts" in name)
        and ("Open MPI" in args.mpi_version or "Intel" in args.mpi_version)
    ):
        print(f"# Skipped due to MPI version\n", file=args.output)
        continue

    dir = test_dir(args.run_dir, name)

    dir.mkdir(parents=True, exist_ok=True)

    labels = get_labels(name)
    labels.append(f"size{ranks}")

    testname = f"[=[precice.{name}]=]"

    testargs = []
    if ranks == 1:
        testargs = [
            args.executable.as_posix(),
            f"--run_test={name}",
            "--log_level=message",
        ]
    else:
        if not args.mpi:
            print(f"# Skipped as MPI not available\n", file=args.output)
            continue

        testargs = [args.mpi_exec.as_posix(), args.mpi_np, str(ranks)]

        if args.mpi_extra:
            testargs.append(args.mpi_extra)

        if "Intel" not in args.mpi_version and "Microsoft" not in args.mpi_version:
            testargs.extend(["--map-by", ":OVERSUBSCRIBE"])

        testargs.extend(
            [args.executable.as_posix(), f"--run_test={name}", "--log_level=message"]
        )

    testargs = [f'"{arg}"' for arg in testargs]
    add_command(args.output, "add_test", [testname] + testargs)

    pargs = [
        testname,
        "PROPERTIES",
        "ENVIRONMENT",
        '"OMP_NUM_THREADS=2;OMP_PROC_BIND=false"',
        "TIMEOUT",
        '"30"',
        '"SKIP_REGULAR_EXPRESSION"',
        '"error 1175;failed to communicate with smpd manager on"',  # smpd post error
        "WORKING_DIRECTORY",
        f'"{dir.as_posix()}"',
        "LABELS",
        f'"{";".join(labels)}"',
    ]

    if args.paths:
        pargs += ["ENVIRONMENT_MODIFICATION", f"PATH=path_list_append:{args.paths}"]

    add_command(args.output, "set_tests_properties", pargs)
    print(file=args.output)

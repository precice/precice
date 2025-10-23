#! python3

import json
import pathlib
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("milestone")
milestone = parser.parse_args().milestone

REPO_CMD = ("git", "rev-parse", "--show-toplevel")
PR_CMD = (
    "gh",
    "pr",
    "list",
    "--search",
    f'is:pr is:closed milestone:"{milestone}"',
    "--limit",
    "200",
    "--json",
    "number,title",
)

root = pathlib.Path(
    subprocess.run(REPO_CMD, stdout=subprocess.PIPE, encoding="utf8").stdout.strip()
)
prlist = json.loads(
    subprocess.run(
        PR_CMD, cwd=str(root.absolute()), stdout=subprocess.PIPE
    ).stdout.strip()
)

prs = {pr["number"]: pr["title"] for pr in prlist}

print("Missing")
for num, title in prs.items():
    cl = root / "docs" / "changelog" / f"{num}.md"
    if not cl.exists():
        print(f"{num:>4} {title}")

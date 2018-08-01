#!/bin/env python3.6
"""
Assembles a trace file for the chromium tracing tool using multiple log files and prints the result.

The tool is available through chromium browsers (e.g. Google Chrome) using the url chrome://tracing or by using the standalone.
"""
from __future__ import print_function, with_statement

import csv
import json
import argparse
import sys
from io import open

CATEGORIES = {
    "advance": "Main Loop",
    "initialize": "Initialization",
    "initializeData": "Initialization",
    "feedbackMesh": "Mesh",
    "receive global mesh": "Mesh,Communication,Receive",
    "send global mesh": "Mesh,Communication,Send",
    "gather mesh": "Mesh,Communication,Gather",
}
DEFAULT_CATEGORY = "default"


class StoreDictKeyPair(argparse.Action):
    """
    Helper class used as an action in argparse to store a dictionary of the format
    KEY=VALUE and sets the according namespace attribute
    """
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        self._nargs = nargs
        super(StoreDictKeyPair, self).__init__(
            option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        my_dict = {}
        for kv in values:
            k, v = kv.split("=")
            my_dict[k] = v
        setattr(namespace, self.dest, my_dict)


def check_and_parse_args():
    """
    Handles the parsing of the command line arguments
    """
    parser = argparse.ArgumentParser(description="Assembles a trace file for"
                                     "the chromium tracing tool using multiple"
                                     "log files and prints it.")
    parser.add_argument("logs", action=StoreDictKeyPair,
                        nargs="+", metavar="PARTICIPANT=LOGFILE")
    parser.add_argument("-p", "--pretty",  action="store_true",
                        help="Print the JSON in a pretty format.")
    parser.add_argument("-d", "--default",  default=DEFAULT_CATEGORY, metavar="CATEGORY",
                        help="The default category for unknown events. "
                        "Known events are: " + ", ".join(CATEGORIES.keys()))
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(sys.argv[1:])


def build_process_name_entry(name, pid):
    """
    Builds a dictionary entry representing a metadata event to name a process id.
    """
    return {
            "name": "process_name",
            "ph": "M",
            "pid": pid,
            "tid": 0,
            "args": {"name": name}
        }


def build_thread_name_entry(name, pid, tid):
    """
    Builds a dictionary entry representing a metadata event to name a thread id
    local to a process.
    """
    return {
            "name": "thread_name",
            "ph": "M",
            "pid": pid,
            "tid": tid,
            "args": {"name": name}
        }


def main():
    args = check_and_parse_args()

    # The output will be in the JSONArray format described in the specification
    traces = []
    for pid, (participant, file) in enumerate(args.logs.items()):
        # The pid identifies each participant and is used as process id
        traces.append(build_process_name_entry(participant, pid))
        with open(file, newline="") as csvfile:
            ranks = set()
            log = csv.DictReader(csvfile, delimiter=",")
            for row in log:
                rank = row["Rank"]
                if rank not in ranks:
                    traces.append(build_thread_name_entry(
                        "Rank {}".format(rank), pid, rank))
                    ranks.add(rank)
                # The current log format contains begin and end timestamps of
                # events which corresponds to the specified duration events
                event = {
                    "name": row["Name"],
                    "cat": CATEGORIES.get(row["Name"], args.default),
                    "tid": rank,
                    "pid": pid,
                    "ts": row["Timestamp"],
                }
                event["ph"] = "B" if row["State"] == "1" else "E"
                traces.append(event)
    if args.pretty:
        print(json.dumps(traces, indent=2))
    else:
        print(json.dumps(traces))


if __name__ == "__main__":
    main()

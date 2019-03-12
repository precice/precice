#!/usr/bin/env python3
""" Plots Events.log files and shows Events on a timeline for each rank. """

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import argparse, sys

colormap = matplotlib.cm.get_cmap()

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description = "Visualize Event logs")
parser.add_argument('groupby', choices = ['name', 'rank'], help = "Group by event or rank")
parser.add_argument('--file', help = "File name of log", type = str, default = "Events.log")
parser.add_argument('--filter', help = "Filter expression used on pandas.query", type = str)
parser.add_argument('--runindex', help = "Index of run, -1 is latest", type = int, default = -1)

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

df = pd.read_csv(args.file, index_col = 0, parse_dates = [0])
df.sort_index()

# Get one dataset (last by default)
df = df.loc[df.index.unique()[args.runindex]]

# We do not care about paused state and treat them as stopped
df.State[df.State == 2] = 0

# Zero timestamps on first event
df.Timestamp = df.Timestamp - min(df.Timestamp)

# Remove the _GLOBAL event,
df = df[df.Name != "_GLOBAL"]

# Apply filter
if args.filter:
    try:
        df.query(args.filter, inplace = True)
    except:
        print("There is something wrong with the filter string you supplied.")
        sys.exit(-1)


height = 1
padding = 0.2

y_pos = []
lefts = []
widths = []
heights = []
names = []
colors = []

if args.groupby == "name":
    groups = ["Name", "Rank"]
elif args.groupby == "rank":
    groups = ["Rank", "Name"]

for i, outer_key in enumerate(df.groupby(groups[0]).groups.keys()):
    ev_df = df[df[groups[0]] == outer_key]
    inner_keys = ev_df[groups[1]].unique()
    sub_height = height / len(inner_keys)

    for j, inner_key in enumerate(inner_keys):
        start = 0
        state = 0
        for index, cols in ev_df[ev_df[groups[1]] == inner_key].iterrows():
            if cols.State == 1 and state == 0: # Event goes from stopped to started
                start = cols.Timestamp
                state = 1
            if cols.State == 0 and state == 1: # Event goes from started to stopped
                lefts.append(start)
                widths.append(cols.Timestamp - start)
                heights.append(sub_height)
                y_pos.append(i * height - height/2 + i * padding + sub_height * j)
                colors.append(colormap(j / len(inner_keys)))
                state = 0
                
        names.append(str(outer_key) + " (" + str(inner_key) + ")")


ax = plt.gca()
ax.barh(y_pos, widths, heights, lefts, color = colors, align = 'center')
ax.set_yticks(sorted(list(set(y_pos)))) # Make y_pos unique
ax.set_yticklabels(names)
ax.invert_yaxis()
ax.set_xlabel("Time [ms]")
ax.grid(True)
ax.set_title("Run at " + str(df.index[0]))

plt.show()


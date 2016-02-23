#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pdb
from EventTimings import parseEventlog

record = parseEventlog("EventTimings.log")[-1] # We just take the last of timing blocks
timings = record["timings"]

fig, ax = plt.subplots()
ax.set_ylabel('Times')
plt.title("Run at " + str(record["timestamp"]))

width = 0.4

avgs = timings["avg"]
mins = timings["min"]
maxs = timings["max"]

x = range(len(timings))

plt.bar(x, avgs, width, align="center")
plt.bar(x, maxs-mins, width=width*0.4, bottom=mins, align="center", color="r")
plt.xticks([i for i in x], timings["name"], rotation=15)
plt.grid()

plt.show()

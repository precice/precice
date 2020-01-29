# Python script steering gnuplot to repeatedly read table data and update plot.

import argparse
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style


def makeParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',
                        type=str,
                        help='The name of the file to plot')
    parser.add_argument(
        "-i",
        "--interval",
        dest="interval",
        default=10000,
        type=int,
        help="Replot interval in milli seconds (default: 10000)")
    parser.add_argument(
        "-x",
        "--xcolumn",
        dest="xcol",
        default=0,
        type=int,
        help="Index of data file column holding x values (default: 0)")
    parser.add_argument(
        "-y",
        "--ycolumn",
        dest="ycol",
        default=1,
        type=int,
        help="Index of data file column holding y values (default: 1)")
    parser.add_argument(
        "-l",
        "--last",
        dest="last",
        default=-1,
        type=int,
        help="Only last N lines of data values will be plotted")
    parser.add_argument("-n",
                        "--name",
                        dest="name",
                        default="",
                        type=str,
                        help="Sets the name of the plot (default: no name)")
    parser.add_argument(
        "-e",
        "--every",
        dest="every",
        default=1,
        type=int,
        help="Sets the distance of data lines read (every option) (default: 1)"
    )

    return parser


def loadFile(filename):
    graph_data = open(filename, 'r').read()
    lines = list(filter(None, graph_data.split('\n')))
    header = list(filter(None, lines[0].split(' ')))
    splitLines = []
    for line in lines[1:]:
        cols = list(filter(None, line.split(' ')))
        assert (len(cols) == len(header))
        splitLines.append(cols)
    return header, splitLines


def main():
    args = makeParser().parse_args()

    style.use('fivethirtyeight')

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    if (args.name):
        ax1.set_title(args.name)

    def animate(i):
        header, lines = loadFile(args.filename)
        assert (len(header) >= max(args.xcol, args.ycol))

        if (args.every > 1):
            lines = [
                l
                for (i, l) in enumerate(lines)
                if i % args.every == 0
            ]

        xs, ys = [], []
        for line in lines:
            xs.append(float(line[args.xcol]))
            ys.append(float(line[args.ycol]))

        if (args.last and args.last > 0):
            if len(xs) > args.last: del xs[len(xs) - args.last:]
            assert(len(xs) <= args.last)
            if len(ys) > args.last: del ys[len(ys) - args.last:]
            assert(len(ys) <= args.last)

        ax1.clear()
        ax1.set_xlabel(header[args.xcol])
        ax1.set_ylabel(header[args.ycol])
        ax1.plot(xs, ys)

    animation.FuncAnimation(fig, animate, interval=1000)
    animate(0)
    plt.show()


if __name__ == '__main__':
    main()

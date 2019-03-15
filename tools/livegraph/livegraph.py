# Python script steering gnuplot to repeatedly read table data and update plot.

import os
import sys
from optparse import OptionParser
import time

parser = OptionParser("%prog [options] filename")
parser.add_option ("-i", "--interval", dest="interval", default=10000,
                   type="int", help="Replot interval in milli seconds (default: 10000)")
parser.add_option ("-x", "--xcolumn", dest="xcol", default="0", type="string", 
                   help="Index of data file column holding x values (default: 0)")
parser.add_option ("-y", "--ycolumn", dest="ycol", default="1", type="string", 
                   help="Index of data file column holding y values (default: 1)")
parser.add_option ("-l", "--last", dest="last", default=-1, type="int",
                   help="Only last N lines of data values will be plotted")
parser.add_option ("-t", "--terminal", dest="terminal", default="", type="string",
                   help="Sets the terminal type to plot to (default: GNUPlot default)")
parser.add_option ("-n", "--name", dest="name", default="", type="string",
                   help="Sets the name of the plot (default: no name)"  )
parser.add_option ("-e", "--every", dest="every", default=1, type="int",
                   help="Sets the distance of data lines read (every option) (default: 1)")
(options, args) = parser.parse_args()
if len(args) != 1:
    parser.print_help()
    sys.exit()
filename = args[0]
gnuplot = os.popen ( 'gnuplot', 'w' )
gnuplot.write ( 'set title "' + options.name + '";' )
if options.terminal != "":
    gnuplot.write ( 'set terminal ' + options.terminal + ';' )
plotstring = 'plot "'
if options.last != -1:
    plotstring += '< tail -' + str(options.last) + ' '
plotstring += filename + '" every ' + str(options.every) + ' using ' + \
              options.xcol + ':' + options.ycol + ' with lines lw 2\n'
print 'Plot command:', plotstring
gnuplot.write ( plotstring )
gnuplot.write ( 'set autoscale;' )
gnuplot.flush ()
try:
    while True:
        time.sleep ( options.interval / 1000.0 )
        gnuplot.write ( 'replot\n' )
        gnuplot.flush ()
finally:
    gnuplot.close ()

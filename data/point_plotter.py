#!/bin/env python

import getopt, sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

from mpl_toolkits import basemap

##
#  @class Plot
#  @brief Returns a basic plot given a set of parameters.
#
#  This forms the basis for generating a plot using this suite of tools.
#  It returns a Matplotlib plot with certain parameters already set up.
class Plot:
    
    ##
    #  Initializes the Plot with a set of basic parameters. Use the
    #  addsubplot() method to add a sub-plot to the plot.
    #
    #  @param title The title for the plot.
    #  @param xlabel The label to be displayed on the x-axis.
    #  @param ylabel The label for the y-axis.
    #  @param legend A legend to be displayed on the lower left of the plot.
    #  @param width The width of the plot in inches (dpi = 100).
    #  @param height The height of the plot in inches (dpi = 100).
    def __init__(self, title = None, xlabel = None, ylabel = None, legend = None, width = 10, height = 10):
        ## Defines the figure object to which we can add subplots.
        self.figure = plt.figure(figsize=(width, height), dpi=100)
        
        if ylabel != None:
            plt.ylabel(ylabel, fontsize=14)
        
        if xlabel != None:
            plt.xlabel(xlabel, fontsize=14)
        
        if title != None:
            plt.title(title)

        if legend != None:
            plt.legend(legend, loc='lower left')

        ## Internal counter for how many subplots we have.
        self.subplotcounter = 1

    ##
    #  Adds a subplot to the figure and returns it.
    #
    #  @return The subplot that has been added to the already generated plot.
    def addsubplot(self):
        retval = self.figure.add_subplot(1, 1, 1)
        self.subplotcounter += 1;
        return retval
    
    ##
    #  Shows the plot.
    def show(self):
        plt.show()
    
    ##
    #  Saves the figure to disk.
    #
    #  @param filename The name fo the file to save.
    def savefig(self, filename):
        plt.savefig(filename)

def main():
    opts, args = getopt.getopt(sys.argv[1:], "f:x:y:", ["file", "x", "y"])
    
    files = None
    x = None
    y = None
    colors = ["r", "b", "g"]
    
    for o, a in opts:
        if o in ("-f", "--file"):
            files = a.split(",")
        elif o in ("-x", "--x"):
            x = a
        elif o in ("-y", "--y"):
            y = a
        elif o in ("-z", "--z"):
            z = a

    if files == None:
        print("No file selected")
        sys.exit(-1)

    data = {}
    plot_index = 0

    for file in files:
        file = file.strip()
        data[plot_index] = {}
        with open(file, "r") as f:
            for line in f:
                split_line = line.split()
                if x == split_line[0] and y == split_line[1]:
                    data[plot_index][101 - int(split_line[2])] = split_line[4]
        plot_index += 1

    p = Plot("CCA Depth Profile at (%s, %s)" % (x, y), "Units", "Depth (Grid Points)", None, 7, 10)

    plot_index = 0
    for file in files:
        xvals = []
        yvals = []
        file = file.replace("CCA0", "")
        file = file.replace(".ascii", "")
        for i in xrange(1, 100):
            xvals.append(data[plot_index][i])
            yvals.append(i)
        p.addsubplot().plot(xvals, yvals, "-", color=colors[plot_index], label="Iteration " + file)
        plot_index += 1

    plt.legend(loc="lower left")
    plt.ylim(plt.ylim()[::-1])
    plt.show()

    return 0

if __name__ == "__main__":
    main()
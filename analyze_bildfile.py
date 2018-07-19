#!/usr/bin/env python
# analyze the bild files form relion to get an idea of how even the representation of views is

import numpy as np
import sys
import matplotlib.pyplot as plt

def calc_lengths(line):
    linesplit = line.split()
    x1 =float(linesplit[1])
    y1 =float(linesplit[2])
    z1 =float(linesplit[3])
    x2 =float(linesplit[4])
    y2 =float(linesplit[5])
    z2 =float(linesplit[6])
    distance = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    return(distance)

def check_bildfile(bildfile):
    file = open(bildfile,'r').readlines()
    lengths = []
    for i in file:
        if '.cylinder' in i:
            lengths.append(calc_lengths(i))
    std,n,min,max = np.std(lengths),len(lengths),np.min(lengths),np.max(lengths)
    return(lengths)


nplots = len(sys.argv[1:])

fig = plt.figure()


# get the y limit in a very clunky way :()
ylims = []
for i in sys.argv[1:]: 
    plt.hist(check_bildfile(i))
    ymin,ymax = plt.ylim()
    ylims.append(ymax)
plt.close()
ymax = max(ylims)

plotcount = 1
for i in sys.argv[1:]:    
    plt.subplot(nplots,1,plotcount)
    plt.title(i)
    plt.ylim(0, ymax)
    plt.ylabel('particles', fontsize=8)
    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.hist(check_bildfile(i))
    plotcount +=1
    plt.tight_layout()
plt.savefig('bildcomp.png')
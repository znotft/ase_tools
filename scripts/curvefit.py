#!/usr/bin/env python

# -*- coding: utf-8 -*-

import numpy as np
import sys

inputfile = sys.argv[1].lstrip()
order = int(sys.argv[2].lstrip())

def read_data(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    xy=[]
    for line in lines:
        xy.append([float(f) for f in line.split()])
    xy = np.array(xy)
    return [xy[:,0], xy[:,1]]
 
def get_parameters(filename, order):
    [x, y] = read_data(filename)
    return np.polyfit(x, y, order, rcond=None, full=False)

def minimum(filename, order):
    [x, y] = read_data(filename)
    p = get_parameters(filename, order)
    ynew = np.polyval(p, np.arange(x.min(), x.max(),0.0002))
    return [np.arange(x.min(), x.max(),0.0002)[ynew.argmin()], ynew.min()]

print '%i degree polynomial:' % order
print 'Minimum at: x: %.6f, y: %.6f' % (minimum(inputfile, order)[0], minimum(inputfile, order)[1] )
print 'Parameters'
p = get_parameters(inputfile, order)
for n in range(order+1):
    print 'a%i: %f' % (n, p[order-n])

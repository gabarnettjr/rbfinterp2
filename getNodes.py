#!/usr/bin/python
"""
Starting from a Cartesian grid, this script will jostle the nodes by some
chosen proportion to create "random" nodes for testing.

Greg Barnett
December 2022
"""
################################################################################

import os
from sys import argv, path

path.append(".")
import IO
import rbf2

################################################################################

# Process the input, if there is any.

coordsDir = "randomCoords"
nx = 16
ny = 16
alp = .3
a = 0
b = 1
c = 0
d = 1

argv = argv[1:]
if len(argv) > 0 :
    coordsDir = argv[0];  argv = argv[1:]
    if not os.path.isdir(coordsDir) :
        s = "First input must be a coordinates directory."
        raise ValueError(s)

if len(argv) > 0 :  nx =   int(argv[0]);  argv = argv[1:]
if len(argv) > 0 :  ny =   int(argv[0]);  argv = argv[1:]
if len(argv) > 0 : alp = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   a = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   b = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   c = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   d = float(argv[0]);  argv = argv[1:]

################################################################################

# Jostle the nodes and save them in the coordinates directory.

xx, yy = rbf2.jostle(nx, ny, alp, a, b, c, d)

IO.saveArray(os.path.join(coordsDir, "x.txt"), xx)
IO.saveArray(os.path.join(coordsDir, "y.txt"), yy)


#!/usr/bin/python
"""
Sometimes you already have some nodes that you are working with, (x, y), and
what you need are some evaluation points (xe, ye) to interpolate to.  This
script will read the files x.txt and y.txt, and use their values to generate
a nice collection of corresponding evaluation points for testing.

Greg Barnett
January 2023
"""
################################################################################

import os
from sys import argv, path
import numpy as np

path.append(".")
import IO
import rbf2

################################################################################

# Process the input, if there is any.

coordsDir = "randomCoords"
nx = 32
ny = 32
alp = 0
a = b = c = d = ""

argv = argv[1:]
if len(argv) > 0 :
	coordsDir = argv[0];  argv = argv[1:]
	if not os.path.isdir(coordsDir) :
		s = "First input must be a coordinates directory."
		raise ValueError(s)
	elif (not os.path.isfile(os.path.join(coordsDir, "x.txt"))) \
    or   (not os.path.isfile(os.path.join(coordsDir, "y.txt"))) :
		s = "Coordinates directory must contain nodes."
		raise ValueError(s)

if len(argv) > 0 :  nx =   int(argv[0]);  argv = argv[1:]
if len(argv) > 0 :  ny =   int(argv[0]);  argv = argv[1:]
if len(argv) > 0 : alp = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   a = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   b = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   c = float(argv[0]);  argv = argv[1:]
if len(argv) > 0 :   d = float(argv[0]);  argv = argv[1:]

################################################################################

if (a == "") or (b == "") or (c == "") or (d == "") :
	# Use nodes to get boundaries of eval pts.
	x = IO.loadArray(coordsDir + "\\x.txt")
	y = IO.loadArray(coordsDir + "\\y.txt")
	a = np.min(x)
	b = np.max(x)
	c = np.min(y)
	d = np.max(y)

################################################################################

# Jostle the nodes and then save them.

xe, ye = rbf2.jostle(nx, ny, alp, a, b, c, d)

IO.saveArray(os.path.join(coordsDir, "xe.txt"), xe)
IO.saveArray(os.path.join(coordsDir, "ye.txt"), ye)


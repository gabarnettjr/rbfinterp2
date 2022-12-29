
from sys import argv
import numpy as np
import IO

if len(argv) == 2 :
    wdir = argv[1]
else :
    wdir = "cdCoords\\smoothData"

################################################################################

x = IO.loadArray(wdir + "\\..\\x.txt")
y = IO.loadArray(wdir + "\\..\\y.txt")

xe = IO.loadArray(wdir + "\\..\\xe.txt")
ye = IO.loadArray(wdir + "\\..\\ye.txt")

################################################################################

def eff(x,y):
    return np.cos(2*np.pi*x/25000) * np.sin(2*np.pi*y/25000)

# def eff(x,y):
    # return np.cos(2*np.pi*x/25000) * np.cos(2*np.pi*y/25000) \
    # + np.sin(2*np.pi*x/25000) * np.sin(2*np.pi*y/25000)

# def eff(x,y):
    # return np.cos(2*np.pi*x/50000) * np.cos(2*np.pi*y/50000) \
    # + np.sin(2*np.pi*x/50000) * np.sin(2*np.pi*y/50000)

################################################################################

IO.saveArray(wdir + "\\f.txt", eff(x, y))

IO.saveArray(wdir + "\\fe.txt", eff(xe, ye))


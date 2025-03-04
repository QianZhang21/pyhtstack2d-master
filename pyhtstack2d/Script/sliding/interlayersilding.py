import numpy as np
import os
import shutil
from pyhtstack2d.buildbilayer.RotMovePOSCAR import move_poscar
from pyhtstack2d.buildbilayer.stackBilayer import Bilayer


def interpolate_positions(inipos, finalpos, numpoints):
    return [inipos + (finalpos - inipos) * i / numpoints for i in range(numpoints+1)]


monofile = "FeCl2-1.vasp"
inipos = np.array([1/3, -1/3, 0])
finalpos = np.array([2/3, -2/3, 0])
numpoints = 9

interpos = interpolate_positions(inipos, finalpos, numpoints)
for pos in interpos:
    if os.path.exists("POSCAR_moved"):
        shutil.rmtree("POSCAR_moved")
    move_poscar(monofile, pos)
    movemonofile = os.path.join("POSCAR_moved", os.listdir("POSCAR_moved")[0])
    Bilayer(movemonofile, monofile, overwrite=False, skip_xy_rev=True).WritePOSCAR()

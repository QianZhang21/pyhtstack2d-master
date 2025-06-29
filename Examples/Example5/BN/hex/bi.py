import os
import shutil
from pyhtstack2d.buildbilayer.batchStackBilayer import GenBiLayer
from pyhtstack2d.buildbilayer.RotMovePOSCAR import rotate_poscar, move_poscar


mono = "BN-POSCAR"
if os.path.exists("POSCAR_moved"):
    # os.removedirs("POSCAR_moved")
    shutil.rmtree("POSCAR_moved")
os.makedirs("POSCAR_moved")
if os.path.exists("POSCAR_rotated"):
    # os.removedirs("POSCAR_rotated")
    shutil.rmtree("POSCAR_rotated")
os.makedirs("POSCAR_rotated")
move_poscar(mono, [0, 0, 0])
move_poscar(mono, [1 / 3, -1 / 3, 0])
move_poscar(mono, [-1 / 3, 1 / 3, 0])
rotate_poscar(mono, 0)
rotate_poscar(mono, 60)
formula_w = GenBiLayer(pos_dir="POSCAR_moved", pos_dir2="POSCAR_rotated", genmode="bilayer",
                       overwrite=False, skip_xy_rev=True).batch_stack()


shutil.rmtree("POSCAR_moved")
shutil.rmtree("POSCAR_rotated")



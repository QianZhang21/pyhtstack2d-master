# Import necessary modules
import os
import shutil
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Import custom tools for bilayer construction
from pyhtstack2d.buildbilayer.stackBilayer import Bilayer
from pyhtstack2d.buildbilayer.RemvDuplicates import remove_duplicates
from pyhtstack2d.buildbilayer.builsupercell import build_supercell
from pyhtstack2d.buildbilayer.CenterZ import center


# Swap lattice a and b if a < b to ensure lattice ordering
def swap_a_b(stfile):
    st = Structure.from_file(stfile)
    la = st.lattice.a
    lb = st.lattice.b
    if la < lb:
        lc = st.lattice.matrix[2][2]
        position = st.frac_coords
        newcell = [[lb, 0, 0], [0, la, 0], [0, 0, lc]]
        newposition = np.ones_like(position)
        newposition[:, 0] = position[:, 1]
        newposition[:, 1] = position[:, 0]
        newposition[:, 2] = position[:, 2]
        newst = Structure(lattice=newcell, species=st.species, coords=newposition, coords_are_cartesian=False)
        newst.to(stfile, fmt="poscar")  # Save updated structure
        return newst.lattice.a, newst.lattice.b
    else:
        return la, lb


# Calculate lattice mismatch ratio
def calmismatch(la1, la2):
    return 2 * abs(la1 - la2) / (la1 + la2)


# Extract lattice information from a BN structure file
def BNlattice(BN_file):
    BN_st = Structure.from_file(BN_file)
    BN_num = BN_st.num_sites
    la1_tmp = BN_st.lattice.a
    lb1_tmp = BN_st.lattice.b
    la1 = la1_tmp if la1_tmp > lb1_tmp else lb1_tmp
    lb1 = lb1_tmp if la1_tmp > lb1_tmp else la1_tmp
    return {"lalb": [la1, lb1], "BN_file": BN_file, "BN_num": BN_num}


# Define allowed BN structures with precomputed supercells
BN_allowla = {}
BN_allowla.update({"1sqrt3": BNlattice(BN_file="BN/orh/BN-sqrt3-1.vasp")})
BN_allowla.update({"sqrt7sqrt21": BNlattice(BN_file="BN/orh/BN-sqrt21-sqrt7.vasp")})
BN_allowla.update({"sqrt13sqrt39": BNlattice(BN_file="BN/orh/BN-sqrt13-sqrt39.vasp")})

count = 0
lattice_system = ["orthorhombic", "tetragonal", "cubic"]  # Target lattice systems
dir_path = "surface_dir_square"  # or "surface_dir_square_layer2"
stfilelist = []

# Traverse metal surface structure files
for filepath in os.walk(dir_path):
    for file in filepath[2]:
        stfilelist.append(file)
        st = Structure.from_file(os.path.join(dir_path, file))
        spg_analyzer_st = SpacegroupAnalyzer(st, symprec=0.5, angle_tolerance=5)
        if spg_analyzer_st.get_lattice_type() in lattice_system:
            metal = Structure.from_file(os.path.join(dir_path, file))
            metal_num = metal.num_sites
            la2 = metal.lattice.a
            lb2 = metal.lattice.b

            # Normalize a/b order for consistency
            if la2 < lb2:
                la2, lb2 = swap_a_b(os.path.join(dir_path, file))

            gamma = metal.lattice.gamma
            if abs(gamma - 90) > 1:
                print(f"gamma is not 90, {gamma}")
                continue

            # Try matching with all allowed BN lattices
            super_combine = {}
            for BN, BNinfo in BN_allowla.items():
                BN_file = BNinfo["BN_file"]
                la1 = BNinfo["lalb"][0]
                lb1 = BNinfo["lalb"][1]
                super_combine[BN] = {'la': [], 'lb': [], 'la_mismatch': [], 'lb_mismatch': []}

                # Check matching for a-axis
                for i in range(1, 10):
                    for j in range(1, 10):
                        mismatch_la = calmismatch(la1 * i, la2 * j)
                        if mismatch_la < 0.03:
                            super_combine[BN]['la'].append([i, j])
                            super_combine[BN]['la_mismatch'].append(mismatch_la)

                # Check matching for b-axis
                for i in range(1, 10):
                    for j in range(1, 10):
                        mismatch_lb = calmismatch(lb1 * i, lb2 * j)
                        if mismatch_lb < 0.03:
                            super_combine[BN]['lb'].append([i, j])
                            super_combine[BN]['lb_mismatch'].append(mismatch_lb)

            # Choose combination with minimum number of atoms
            mincombine_BN = None
            mincombine_BN_super = None
            mincombine_metal_super = None
            mincombine_numatoms = 100000
            la_mismatch = None
            lb_mismatch = None

            for sc, scinfo in super_combine.items():
                la = scinfo['la']
                lb = scinfo['lb']
                if len(la) == 0 or len(lb) == 0:
                    continue
                BN_super = la[0][0] * lb[0][0]
                metal_super = la[0][1] * lb[0][1]
                numatoms = BN_allowla[sc]["BN_num"] * BN_super + metal_num * metal_super
                if numatoms < mincombine_numatoms:
                    mincombine_numatoms = numatoms
                    mincombine_BN = sc
                    mincombine_BN_super = [la[0][0], lb[0][0], 1]
                    mincombine_metal_super = [la[0][1], lb[0][1], 1]
                    la_mismatch = scinfo['la_mismatch'][0]
                    lb_mismatch = scinfo['lb_mismatch'][0]

            # If a match is found, build the heterostructure
            if mincombine_BN is not None and mincombine_numatoms < 200:
                count += 1
                print(f"File: {file}, BN: {mincombine_BN}, Supercell: {mincombine_BN_super}, "
                      f"Metal Supercell: {mincombine_metal_super}, la_mismatch: {la_mismatch * 100}, "
                      f"lb_mismatch: {lb_mismatch * 100}, Total Atoms: {mincombine_numatoms}")

                # Build supercells for BN and metal
                if os.path.exists("supercell"):
                    shutil.rmtree("supercell")
                os.mkdir("supercell")
                build_supercell(BN_allowla[mincombine_BN]["BN_file"], mincombine_BN_super, save_path="supercell/BN-POSCAR")
                build_supercell(os.path.join(dir_path, file), mincombine_metal_super, save_path="supercell/metal-POSCAR")

                # Construct bilayer heterostructure
                bigen = Bilayer("supercell/metal-POSCAR", "supercell/BN-POSCAR",
                                overwrite=False, skip_xy_rev=True, d_inter=3.0, lv=28)
                dirpath_w = bigen.formula_w
                ex_count = 0
                while os.path.exists(os.path.join("BiPOSCAR_dir", dirpath_w)):
                    ex_count += 1
                    bigen.formula_w = bigen.formula_w + f"_{ex_count}"
                    dirpath_w = bigen.formula_w
                bigen.WritePOSCAR()

                # Remove duplicates and post-process generated files
                remove_duplicates(os.path.join("BiPOSCAR_dir", dirpath_w))
                for dirpath, dirnames, filenames in os.walk(os.path.join("BiPOSCAR_dir", dirpath_w)):
                    for filename in filenames:
                        if filename.endswith("POSCAR"):
                            stfile = os.path.join(dirpath, filename)
                            st = Structure.from_file(stfile)
                            Bilayer(stfile, "supercell/BN-POSCAR", overwrite=False, savepath=dirpath,
                                    skip_xy_rev=True, savenamemode=2, d_inter=3.0, lv=28).WritePOSCAR()
                            remove_duplicates(dirpath, multist=False)
                            os.remove(stfile)

                            rmpos = []
                            for subdirpath, subdirnames, subfilenames in os.walk(dirpath):
                                for subfilename in subfilenames:
                                    stfile = os.path.join(subdirpath, subfilename)
                                    st = Structure.from_file(stfile)
                                    positions = st.cart_coords
                                    elements = [el.symbol for el in st.species]
                                    sortindex = np.argsort(positions[:, 2])
                                    sorted_elem = [elements[i] for i in sortindex]

                                    # Check if B and N atoms are on both sides of the bilayer (top and bottom)
                                    if sorted_elem[0] in ["B", "N"] and sorted_elem[-1] in ["B", "N"]:
                                        os.rename(stfile, os.path.join(subdirpath, "POSCAR"))
                                        center(os.path.join(subdirpath, "POSCAR"))
                                        print(os.path.join(subdirpath, "POSCAR"),
                                              max(positions[:, 2]) - min(positions[:, 2]))
                                    else:
                                        rmpos.append(stfile)
                            for rmfile in rmpos:
                                os.remove(rmfile)

# Print final count of successful bilayer structures
print(count)

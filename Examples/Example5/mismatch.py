# Import necessary modules
import os
import shutil
import numpy as np
from pymatgen.core import Element
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Import custom modules from pyhtstack2d for bilayer operations
from pyhtstack2d.buildbilayer.stackBilayer import Bilayer
from pyhtstack2d.buildbilayer.RemvDuplicates import remove_duplicates
from pyhtstack2d.buildbilayer.builsupercell import build_supercell
from pyhtstack2d.buildbilayer.CenterZ import center

# Define scaling matrices for 120° and 60° lattice gamma angles
scaling_matrix_dict_120 = {
    "1": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "sqrt3": [[1, -1, 0], [1, 2, 0], [0, 0, 1]],
    "2": [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
    "sqrt7": [[2, -1, 0], [1, 3, 0], [0, 0, 1]],
    "3": [[3, 0, 0], [0, 3, 0], [0, 0, 1]],
    "sqrt12": [[2, -2, 0], [2, 4, 0], [0, 0, 1]],
    "sqrt13": [[3, -1, 0], [1, 4, 0], [0, 0, 1]],
    "4": [[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    "sqrt21": [[4, -1, 0], [1, 5, 0], [0, 0, 1]],
    "5": [[5, 0, 0], [0, 5, 0], [0, 0, 1]],
    "sqrt27": [[3, -3, 0], [3, 6, 0], [0, 0, 1]],
    "sqrt31": [[5, -1, 0], [1, 6, 0], [0, 0, 1]]
}

scaling_matrix_dict_60 = {
    "1": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    "sqrt3": [[1, 1, 0], [-1, 2, 0], [0, 0, 1]],
    "2": [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
    "sqrt7": [[2, 1, 0], [-1, 3, 0], [0, 0, 1]],
    "3": [[3, 0, 0], [0, 3, 0], [0, 0, 1]],
    "sqrt12": [[2, 2, 0], [-2, 4, 0], [0, 0, 1]],
    "sqrt13": [[3, 1, 0], [-1, 4, 0], [0, 0, 1]],
    "4": [[4, 0, 0], [0, 4, 0], [0, 0, 1]],
    "sqrt21": [[4, 1, 0], [-1, 5, 0], [0, 0, 1]],
    "5": [[5, 0, 0], [0, 5, 0], [0, 0, 1]],
    "sqrt27": [[3, 3, 0], [-3, 6, 0], [0, 0, 1]],
    "sqrt31": [[5, 1, 0], [-1, 6, 0], [0, 0, 1]]
}

# Define scale factors for supercell matching
allowsuper = {
    "1": 1, "sqrt3": np.sqrt(3), "2": 2, "sqrt7": np.sqrt(7), "3": 3, "sqrt12": np.sqrt(12),
    "sqrt13": np.sqrt(13), "4": 4, "sqrt21": np.sqrt(21), "5": 5, "sqrt27": np.sqrt(27),
    "sqrt31": np.sqrt(31)
}

# Function to calculate lattice mismatch
def calmismatch(la1, la2):
    return 2 * abs(la1 - la2) / (la1 + la2)

# Load BN structure and get its lattice constant
BN_file = "BN/hex/BN-POSCAR-gamma60"
BN = Structure.from_file(BN_file)
la1 = BN.lattice.a
BN_num = BN.num_sites

# Directory containing surface metal structures
dir_path = "surface_dir"
stfilelist = []
count = 0

# Iterate through all structure files in the directory
for filepath in os.walk(dir_path):
    for file in filepath[2]:
        stfilelist.append(file)
        st = Structure.from_file(os.path.join(dir_path, file))

        # Analyze symmetry and lattice
        spg_analyzer_st = SpacegroupAnalyzer(st, symprec=0.5, angle_tolerance=5)
        spg_st = spg_analyzer_st.get_space_group_number()

        # Process only hexagonal structures
        if spg_analyzer_st.get_lattice_type() == "hexagonal":
            metal = Structure.from_file(os.path.join(dir_path, file))
            metal_num = metal.num_sites
            la2 = metal.lattice.a
            gamma = metal.lattice.gamma

            # Choose appropriate scaling matrix based on gamma angle
            if abs(gamma - 120) < 1.0:
                BN_file = "BN/hex/BN-POSCAR"
                scaling_matrix_dict = scaling_matrix_dict_120
            elif abs(gamma - 60) < 1.0:
                BN_file = "BN/hex/BN-POSCAR-gamma60"
                scaling_matrix_dict = scaling_matrix_dict_60
            else:
                print(f"gamma is not 60 or 120, {gamma}")
                continue

            # Try to find matching supercells
            l1_super = None
            l2_super = None
            super_combine = []
            mincombine_numatoms = 100000

            for i, s1 in allowsuper.items():
                for j, s2 in allowsuper.items():
                    delta_la = calmismatch(la1 * s1, la2 * s2)
                    if delta_la < 0.03:  # Match within 3%
                        super_combine.append([i, j])
                        total_num = BN_num * s1 * s1 + metal_num * s2 * s2
                        if total_num < mincombine_numatoms:
                            l1_super = i
                            l2_super = j
                            mincombine_numatoms = total_num

            if l1_super is None or mincombine_numatoms > 200:
                print(f"no suitable supercell for {file}")
                continue

            count += 1

            # Prepare supercell directory
            if os.path.exists("supercell"):
                shutil.rmtree("supercell")
            os.mkdir("supercell")

            # Build BN and metal supercells
            build_supercell(BN_file, scaling_matrix_dict[l1_super], save_path="supercell/BN-POSCAR")
            build_supercell(os.path.join(dir_path, file), scaling_matrix_dict[l2_super], save_path="supercell/metal-POSCAR")

            # Create bilayer structure
            bigen = Bilayer("supercell/metal-POSCAR", "supercell/BN-POSCAR", overwrite=False, skip_xy_rev=True, d_inter=3.0, lv=25)
            dirpath_w = bigen.formula_w
            ex_count = 0
            while os.path.exists(os.path.join("BiPOSCAR_dir", dirpath_w)):
                ex_count += 1
                bigen.formula_w = bigen.formula_w + f"_{ex_count}"
                dirpath_w = bigen.formula_w
            bigen.WritePOSCAR()

            # Remove duplicate structures
            remove_duplicates(os.path.join("BiPOSCAR_dir", dirpath_w))

            # Further clean and center the bilayer structures
            for dirpath, dirnames, filenames in os.walk(os.path.join("BiPOSCAR_dir", dirpath_w)):
                for filename in filenames:
                    if filename.endswith("POSCAR"):
                        stfile = os.path.join(dirpath, filename)
                        Bilayer(stfile, "supercell/BN-POSCAR", overwrite=False, savepath=dirpath,
                                skip_xy_rev=True, savenamemode=2, d_inter=3.0, lv=25).WritePOSCAR()
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

                                # Check if B and N are at top and bottom (indicating correct BN sandwiching)
                                if sorted_elem[0] in ["B", "N"] and sorted_elem[-1] in ["B", "N"]:
                                    os.rename(stfile, os.path.join(subdirpath, "POSCAR"))
                                    center(os.path.join(subdirpath, "POSCAR"))
                                    print(os.path.join(subdirpath, "POSCAR"))
                                else:
                                    rmpos.append(stfile)

                        # Remove invalid structures
                        for rmfile in rmpos:
                            os.remove(rmfile)

# Print total number of valid bilayer structures created
print(count)

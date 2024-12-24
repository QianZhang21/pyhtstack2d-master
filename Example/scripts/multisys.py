# =============================================================================
from pyhtstack2d.tools.queryDB import getPOSACR, getuid

db = "../db_get/c2db-2022-11-30.db"
criteria = 'gap>0, natoms<=10, magstate=NM, thermodynamic_stability_level=3, ' \
           'dynamic_stability_phonons=high, dynamic_stability_stiffness=high'
# uid = getuid(db, criteria, save_csv=True)
# print(len(uid))
getPOSACR(db, criteria, overwrite=False)
# =============================================================================


# =============================================================================
import os
from pyhtstack2d.buildbilayer.stackBilayer import Monolayer

posdir = "POSCAR_dir_c2db"

for poscar in os.listdir(posdir):
    st = Monolayer(os.path.join(posdir, poscar))
    cos_angle1 = abs(st.cos_angle1)
    if abs(cos_angle1 - 0.5) > 0.01 and abs(cos_angle1) > 0.01:
        print(poscar, st.a, st.b, st.cos_angle1, st.cos_angle2, st.cos_angle3)
        os.remove(os.path.join(posdir, poscar))
# =============================================================================


# =============================================================================
from pyhtstack2d.tools.genInput import GenRunDir
GenRunDir(genmode="vaspkit").genInputfile()
# =============================================================================


# =============================================================================
import json
import os
from pymatgen.core import Structure

def swap_a_b_in_structure(structure):
    """
    Swap the lattice vectors a and b in a pymatgen Structure object and adjust coordinates accordingly.

    Parameters:
        structure (Structure): A pymatgen Structure object.

    Returns:
        Structure: A new Structure object with swapped lattice vectors and adjusted coordinates.
    """
    # Get the lattice vectors
    lattice = structure.lattice
    a, b, c = lattice.matrix

    # Swap lattice vectors a and b
    new_lattice = [b, a, c]

    # Adjust fractional coordinates by swapping the first two columns
    new_coords = []
    for frac_coords in structure.frac_coords:
        new_coords.append([frac_coords[1], frac_coords[0], frac_coords[2]])

    # Create a new structure with the updated lattice and coordinates
    new_structure = Structure(new_lattice, structure.species, new_coords, coords_are_cartesian=False)

    return new_structure

with open("monoinfo.json", "r") as f:
    monoinfo = json.load(f)

for key, item in monoinfo.items():
    st = Structure(item['lattice'], item['elements'], item['cart_coords'], coords_are_cartesian=True)
    if st.lattice.a > st.lattice.b and st.lattice.a - st.lattice.b > 0.01:
        st = swap_a_b_in_structure(st)
    st.to(filename=os.path.join("POSCAR_dir", key.split("/")[0] + "-POSCAR"), fmt="poscar")
# =============================================================================


# =============================================================================
import os
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import shutil

stfilelist = []
for filepath in os.walk("POSCAR_dir"):
    for file in filepath[2]:
        stfilelist.append(file)

systemlist = []
for stfile in stfilelist:
    st = Structure.from_file(os.path.join("POSCAR_dir", stfile))
    spg_analyzer_st = SpacegroupAnalyzer(st, symprec=0.1, angle_tolerance=5)
    spg_st = spg_analyzer_st.get_space_group_number()
    if spg_analyzer_st.get_lattice_type() == "monoclinic" and abs(st.lattice.a - st.lattice.b) < 0.01:
        print(stfile, st.lattice.abc, st.lattice.angles)
        # shutil.move(os.path.join("POSCAR_dir", stfile), os.path.join("Orthorhombic_POSCAR_dir", stfile))
    # print(stfile, spg_st, spg_analyzer_st.get_crystal_system())
    systemlist.append(spg_analyzer_st.get_lattice_type())
print(set(systemlist))
# =============================================================================


# =============================================================================
import os
from pymatgen.core import Structure
import numpy as np

# Helper function to move structure
def move_structure_to_origin(st):
    pos = st.frac_coords
    # Calculate the distances to (0, 0, z) for all atoms
    distances = np.linalg.norm(pos[:, :2], axis=1)  # Ignore z-axis, only consider x and y
    distances2 = np.linalg.norm(pos[:, :2]-[1/3, 2/3], axis=1)
    distances3 = np.linalg.norm(pos[:, :2]-[2/3, 1/3], axis=1)

    dmin = min(distances)
    dmin2 = min(distances2)
    dmin3 = min(distances3)

    if dmin < dmin2 and dmin < dmin3:
        closest_atom_idx = np.argmin(distances)  # Find the atom closest to the z-axis
        closest_atom_coords = pos[closest_atom_idx]  # Get the fractional coordinates of the closest atom
        translation_vector = np.array([0, 0, closest_atom_coords[2]]) - closest_atom_coords  # Calculate the translation vector to move the closest atom to (0, 0, z)
    elif dmin2 < dmin and dmin2 < dmin3:
        closest_atom_idx = np.argmin(distances2)
        closest_atom_coords = pos[closest_atom_idx]
        translation_vector = np.array([1/3, 2/3, closest_atom_coords[2]]) - closest_atom_coords
    else:
        closest_atom_idx = np.argmin(distances3)
        closest_atom_coords = pos[closest_atom_idx]
        translation_vector = np.array([2/3, 1/3, closest_atom_coords[2]]) - closest_atom_coords

    # Apply the translation vector to all fractional coordinates
    new_pos = pos + translation_vector

    # Ensure fractional coordinates stay within [0, 1) using modulo operation
    new_pos %= 1

    # Create a new structure with updated positions
    new_structure = Structure(lattice=st.lattice, species=st.species, coords=new_pos, coords_are_cartesian=False)
    return new_structure

# Iterate over files and process structures
full_path_list = []
for filepath in os.walk("Hexagonal_POSCAR_dir"):
    for file in filepath[2]:
        full_path_list.append(os.path.join(filepath[0], file))

for full_path in full_path_list:
    try:
        # Load structure
        st = Structure.from_file(full_path)

        # Get fractional coordinates
        pos = st.frac_coords

        # Check if any atom is at (0, 0, z), (1/3, 2/3, z), or (2/3, 1/3, z) with a tolerance of 0.01
        tolerance = 0.01
        condition1 = np.any(np.all(np.isclose(pos[:, :2], [0, 0], atol=tolerance), axis=1))
        condition2 = np.any(np.all(np.isclose(pos[:, :2], [1 / 3, 2 / 3], atol=tolerance), axis=1))
        condition3 = np.any(np.all(np.isclose(pos[:, :2], [2 / 3, 1 / 3], atol=tolerance), axis=1))
        condition4 = np.any(np.all(np.isclose(pos[:, :2], [0, 1], atol=tolerance), axis=1))
        condition5 = np.any(np.all(np.isclose(pos[:, :2], [1, 0], atol=tolerance), axis=1))
        condition6 = np.any(np.all(np.isclose(pos[:, :2], [1, 1], atol=tolerance), axis=1))

        if not (condition1 or condition2 or condition3):
            # No atom at (0, 0, z), (1/3, 2/3, z), or (2/3, 1/3, z), move structure
            st = move_structure_to_origin(st)
            print("Processing: ", full_path)
            # Save the modified structure back to a file
            os.remove(full_path)
            st.to(fmt="poscar", filename=full_path)
    except Exception as e:
        print(f"Error processing {full_path}: {e}")
# =============================================================================


# =============================================================================
import os
import shutil
import numpy as np
from pyhtstack2d.buildbilayer.batchStackBilayer import GenBiLayer
from pyhtstack2d.buildbilayer.RotMovePOSCAR import rotate_poscar, move_poscar
from pyhtstack2d.buildbilayer.stackBilayer import Monolayer, Bilayer


class Bilayer_modify(Bilayer):
    def check_lattice_mismatch(self, la1, lb1, angle1, la2, lb2, angle2):
        mislattice_a = 2 * abs(la1 - la2) / (la1 + la2) * 100
        mislattice_b = 2 * abs(lb1 - lb2) / (lb1 + lb2) * 100
        self.mismatch = max(mislattice_a, mislattice_b)
        angle1 = -0.5
        la_mean = (la1 + la2) / 2
        lb_mean = la_mean
        # lb_mean = (lb1 + lb2) / 2
        return np.array(
            [[la_mean, 0, 0], [lb_mean * angle1, lb_mean * np.sqrt(1 - angle1 ** 2), 0], [0, 0, self.lv]])

class GenBiLayer_modify(GenBiLayer):
    def get_pos_inf(self, pos_dir):
        pos_obj = []
        la = []
        lb = []
        self.posfilename = []
        if isinstance(pos_dir, str):
            for pos_file in os.listdir(pos_dir):
                self.posfilename.append(pos_file)
                momolayer = Monolayer(os.path.join(pos_dir, pos_file))
                pos_obj.append(momolayer)
                la.append(momolayer.a)
                lb.append(momolayer.b)
        elif isinstance(pos_dir, list):
            for pos_file in pos_dir:
                self.posfilename.append(pos_file)
                momolayer = Monolayer(pos_file)
                pos_obj.append(momolayer)
                la.append(momolayer.a)
                lb.append(momolayer.b)
        else:
            raise ValueError("pos_obj should be a string or a list of strings")
        return pos_obj, la, lb

    def batch_stack(self):
        for mono1_index in self.match_pos_dict.keys():
            mono1 = self.pos_obj[mono1_index]
            if self.single_dir:
                for mono2_index in self.match_pos_dict[mono1_index]:
                    mono2 = self.pos_obj[mono2_index]
                    dz = mono1.get_intra_distances().sum()
                    dz2 = mono2.get_intra_distances().sum()
                    if 30 < dz + dz2 + 15 < 40:
                        lv = 40.0
                    elif 40 < dz + dz2 + 15 < 50:
                        lv = 50.0
                    else:
                        lv = 30.0
                    Bilayer_modify(st1=mono1, st2=mono2, lv=lv, **self.kwargs).WritePOSCAR()
            else:
                for mono2_index in self.match_pos_dict[mono1_index]:
                    mono2 = self.pos_obj2[mono2_index]
                    dz = mono1.get_intra_distances().sum()
                    dz2 = mono2.get_intra_distances().sum()
                    if 30 < dz + dz2 + 15 < 40:
                        lv = 40.0
                    elif 40 < dz + dz2 + 15 < 50:
                        lv = 50.0
                    else:
                        lv = 30.0
                    Bilayer_modify(st1=mono1, st2=mono2, lv=lv, **self.kwargs).WritePOSCAR()

        savedir = Bilayer_modify(self.pos_obj[0], **self.kwargs).savepath
        print(f"Stacking of bilayers is completed. The POSCAR files are saved in the {savedir} directory.")


bi = GenBiLayer_modify("Hexagonal_POSCAR_dir", la_mismatch=0.5)
bidict = bi.match_pos_dict
biobjec = bi.posfilename
for mono1_index in bidict.keys():
    mono1 = os.path.join("Hexagonal_POSCAR_dir", biobjec[mono1_index])
    for mono2_index in bidict[mono1_index]:
        mono2 = os.path.join("Hexagonal_POSCAR_dir", biobjec[mono2_index])
        if os.path.exists("POSCAR_moved"):
            # os.removedirs("POSCAR_moved")
            shutil.rmtree("POSCAR_moved")
        os.makedirs("POSCAR_moved")
        if os.path.exists("POSCAR_rotated"):
            # os.removedirs("POSCAR_rotated")
            shutil.rmtree("POSCAR_rotated")
        os.makedirs("POSCAR_rotated")
        move_poscar(mono1, [0, 0, 0])
        move_poscar(mono1, [1 / 3, -1 / 3, 0])
        move_poscar(mono1, [-1 / 3, 1 / 3, 0])
        rotate_poscar(mono2, 0)
        rotate_poscar(mono2, 60)
        GenBiLayer_modify(pos_dir="POSCAR_moved", pos_dir2="POSCAR_rotated", genmode="bilayer",
                          overwrite=False, skip_xy_rev=True).batch_stack()
# =============================================================================



# =============================================================================
from pyhtstack2d.tools.genInput import GenRunDir
GenRunDir(posdir="BiPOSCAR_dir", multilevel=4, genmode="vaspkit").genInputfile()
# =============================================================================




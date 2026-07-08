import os
import numpy as np
from ase import io
from ase.build import surface
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pyhtstack2d.buildbilayer.RemvDuplicates import remove_duplicates


min_distance_allow = 4.0  # Allowed minimum distance between atoms in the surface structure


def getsuface(stfile, miller_indices, layers, savepath=None, vacuum_layer=5.0):
    structure = io.read(stfile)
    surf = surface(structure, miller_indices, layers, vacuum=vacuum_layer)
    cell = surf.cell
    cell[2, 2] += vacuum_layer
    surf.set_cell(cell, scale_atoms=True)
    mid = "".join(stfile.split('\\')[-1].split('-')[:-1])
    filename = mid + "(" + "".join([str(i) for i in miller_indices]) + ").vasp"
    if savepath:
        filename = os.path.join(savepath, filename)
    surf.write(filename)


def calculate_min_distance(structure):
    # Get the minimum distance between atoms in the structure
    min_distance = float('inf')  # Initialize to infinity

    # Iterate through all pairs of sites in the structure
    for site1 in structure.sites:
        for site2 in structure.sites:
            if site1 != site2:
                # Calculate the distance between the two sites
                distance = site1.distance(site2)
                # Update the minimum distance if this distance is smaller
                if distance < min_distance:
                    min_distance = distance
    if min_distance == float('inf'):
        return min(structure.lattice.a, structure.lattice.b, structure.lattice.c)
    else:
        return min_distance


def compare_distances(surface_dir):
    rm_files = []
    for dir_path, _, files in os.walk(surface_dir):
        for file in files:
            # Check if the file is a structure file
            if file.endswith(('.cif', '.vasp', '.poscar')):
                stfile = os.path.join(dir_path, file)
                # Load the structure from the file
                st = Structure.from_file(stfile)

                # Calculate the minimum distance between atoms in the structure
                min_distance = calculate_min_distance(st)

                if min_distance >= min_distance_allow:
                    bulk = Structure.from_file(os.path.join("3D_POSCAR_dir", file.split("(")[0] + "-POSCAR"))
                    bulk_min_distance = calculate_min_distance(bulk)
                    # Check if the minimum distance is greater than the bulk bond length
                    if min_distance > bulk_min_distance:
                        rm_files.append(stfile)
                        print(f"Remove {stfile} due to minimum distance {min_distance} > {bulk_min_distance}")
                    elif st.num_sites == 1 and (st.lattice.a > bulk_min_distance or st.lattice.b > bulk_min_distance):
                        rm_files.append(stfile)
                        print(
                            f"Remove {stfile} due to single atom with large lattice parameters {st.lattice.a}, {st.lattice.b}")

    for file in rm_files:
        os.remove(file)
        print(f"Remove {file}")


# # ======================================================================================================
surface_dir = 'surface_dir'
layers = 1

if not os.path.exists(surface_dir):
    os.makedirs(surface_dir)

for dir_path, _, files in os.walk('3D_POSCAR_dir'):
    for file in files:
        structure = os.path.join(dir_path, file)
        for miller_indices in [[1, 0, 0], [1, 1, 0], [1, 1, 1],
                               [0, 1, 0], [1, 0, 1],
                               [0, 0, 1], [0, 1, 1]]:
            getsuface(structure, miller_indices, layers=layers, savepath=surface_dir)

# remove_duplicates(surface_dir, multist=False)  # Remove duplicates after generating surfaces

rm_files = []
for dir_path, _, files in os.walk(surface_dir):
    for file in files:
        stfile = os.path.join(dir_path, file)
        st = Structure.from_file(stfile)
        spg_analyzer_st = SpacegroupAnalyzer(st, symprec=0.5, angle_tolerance=5)
        spg_st = spg_analyzer_st.get_space_group_number()
        pos_z = st.cart_coords[:, 2]
        max_z = np.max(pos_z)
        min_z = np.min(pos_z)
        # Check if the lattice type is hexagonal and the thickness is within the allowed range
        if spg_analyzer_st.get_lattice_type() != "hexagonal" or max_z - min_z > 3.2:
            rm_files.append(stfile)
            print(stfile)

for file in rm_files:
    os.remove(file)
    print(f"Remove {file}")

remove_duplicates(surface_dir, multist=False)  # Remove duplicates after filtering surfaces

compare_distances(surface_dir)  # Remove files with minimum distance greater than the bulk bond length
# # ======================================================================================================

# # ======================================================================================================
surface_dir = 'surface_dir_square'
layers = 1
# # For square surface layer 2, change the parameters accordingly
# surface_dir = 'surface_dir_square_layer2'
# layers = 2

if not os.path.exists(surface_dir):
    os.makedirs(surface_dir)

for dir_path, _, files in os.walk('3D_POSCAR_dir'):
    for file in files:
        structure = os.path.join(dir_path, file)
        for miller_indices in [[1, 0, 0], [1, 1, 0], [1, 1, 1],
                               [0, 1, 0], [1, 0, 1],
                               [0, 0, 1], [0, 1, 1]]:
            getsuface(structure, miller_indices, layers=layers, savepath=surface_dir)

lattice_system = ["orthorhombic", "tetragonal", "cubic"]
rm_files = []
for dir_path, _, files in os.walk(surface_dir):
    for file in files:
        stfile = os.path.join(dir_path, file)
        st = Structure.from_file(stfile)
        spg_analyzer_st = SpacegroupAnalyzer(st, symprec=0.5, angle_tolerance=5)
        spg_st = spg_analyzer_st.get_space_group_number()
        pos_z = st.cart_coords[:, 2]
        max_z = np.max(pos_z)
        min_z = np.min(pos_z)
        if spg_analyzer_st.get_lattice_type() not in lattice_system or abs(st.lattice.gamma - 90) > 5 or max_z - min_z > 3.5:
            rm_files.append(stfile)

for file in rm_files:
    os.remove(file)
    print(f"Remove {file}")

remove_duplicates(surface_dir, multist=False)
compare_distances(surface_dir)
# # ======================================================================================================

if os.path.exists("surface_dir_layer2"):
    for dir_path, _, files in os.walk("surface_dir_layer2"):
        for file in files:
            stfile = os.path.join(dir_path, file)
            st = Structure.from_file(stfile)
            min_distance = calculate_min_distance(st)
            if os.path.exists(os.path.join("surface_dir", file)):
                min_distance_single = calculate_min_distance(Structure.from_file(os.path.join("surface_dir", file)))
                if abs(min_distance - min_distance_single) < 0.01:
                    print(f"Remove {stfile}")
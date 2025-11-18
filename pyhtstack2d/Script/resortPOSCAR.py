"""
Uniformly sort POSCAR elements within isomorphous folders.
"""

import os
from pymatgen.core import Structure


def sort_elements_in_poscar(structure, unique_elements):
    """
    This function sorts the elements and their corresponding atomic positions
    in a fixed order provided by unique_elements.
    """
    # Get all elements and positions
    all_elements = structure.species
    all_positions = structure.frac_coords

    # Initialize sorted lists
    sorted_elements = []
    sorted_positions = []

    # For each element in the fixed order (unique_elements),
    # sort positions accordingly
    for element in unique_elements:
        for idx, el in enumerate(all_elements):
            if el == element:
                sorted_elements.append(el)
                sorted_positions.append(all_positions[idx])

    # Create a new structure with the sorted elements and positions
    sorted_structure = Structure(
        lattice=structure.lattice,
        species=sorted_elements,
        coords=sorted_positions,
        coords_are_cartesian=False
    )

    return sorted_structure


def process_poscar_directory(root_dir):
    midlist = os.listdir(root_dir)
    for mid in midlist:
        unique_elements = None
        structures = []  # To store all structures in a directory
        midpath = os.path.join(root_dir, mid)
        stackpatternlist = os.listdir(midpath)

        for stackpattern in stackpatternlist:
            stackpatternpath = os.path.join(midpath, stackpattern)
            cordlist = os.listdir(stackpatternpath)
            for cord in cordlist:
                pospath = os.path.join(stackpatternpath, cord, 'POSCAR')
                structure = Structure.from_file(pospath)
                structures.append(structure)

                # Extract unique elements and sort them by first appearance
                if unique_elements is None:
                    all_elements = []
                    elements = structure.species
                    all_elements.extend(elements)
                    unique_elements = sorted(set(all_elements), key=lambda x: str(x))

                # Sort elements in the current POSCAR structure
                sorted_structure = sort_elements_in_poscar(structure, unique_elements)

                # Save the sorted structure back to POSCAR
                sorted_structure.to(fmt="poscar", filename=os.path.join(stackpatternpath, cord, 'POSCAR'))


# Main program entry, process POSCAR files in the directory
root_directory = 'BiPOSCAR_dir'
process_poscar_directory(root_directory)

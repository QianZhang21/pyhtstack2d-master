import numpy as np
from pymatgen.core import Structure
import warnings


def center(stfile):
    """
    Center the structure in the Z direction.
    """
    # Read the POSCAR file
    structure = Structure.from_file(stfile)
    # Get the lattice vectors
    lattice = structure.lattice
    if lattice.a > lattice.c or lattice.b > lattice.c:
        warnings.warn("The Z direction is not the largest direction.")
    # Get the coordinates of the atoms
    coords = structure.frac_coords
    # Get the Z coordinates of the atoms
    z_coords = coords[:, 2]
    # Get the minimum and maximum Z coordinates
    min_z = np.min(z_coords)
    max_z = np.max(z_coords)
    # Calculate the center of the Z coordinates
    center_z = (min_z + max_z) / 2
    # Calculate the shift needed to center the structure
    shift = 0.5 - center_z
    # Shift the Z coordinates of the atoms
    coords[:, 2] += shift
    # Create a new structure with the shifted coordinates
    new_structure = Structure(lattice=lattice, species=structure.species, coords=coords, coords_are_cartesian=False)
    # Write the new structure to a new POSCAR file
    new_structure.to(stfile, fmt="poscar")

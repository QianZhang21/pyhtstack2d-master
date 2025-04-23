import numpy as np
import os
import shutil
from pyhtstack2d.buildbilayer.RotMovePOSCAR import move_poscar
from pyhtstack2d.buildbilayer.stackBilayer import Bilayer


def interpolate_positions(inipos, finalpos, num_points):
    """
    Generate interpolated positions between the initial and final positions.

    Parameters:
    inipos (numpy.ndarray): Initial position as a numpy array.
    finalpos (numpy.ndarray): Final position as a numpy array.
    num_points (int): Number of interpolation points.

    Returns:
    list: A list of interpolated positions.
    """
    return [inipos + (finalpos - inipos) * i / num_points for i in range(num_points+1)]


# Define initial and final interlayer shift positions
initial_position = np.array([1 / 3, -1 / 3, 0])
final_position = np.array([2 / 3, -2 / 3, 0])

# Number of intermediate sliding steps
num_interpolation_points = 6

# Generate interpolated positions
interpolated_positions = interpolate_positions(initial_position, final_position, num_interpolation_points)


# Define input monolayer structure file
monolayer_file_path = "POSCAR_dir"
moved_dir = "POSCAR_moved"

for monolayer_file in os.listdir(monolayer_file_path):
    monolayer_file = os.path.join(monolayer_file_path, monolayer_file)
    # Loop through each interpolated position and generate bilayer structures
    for position in interpolated_positions:
        # Ensure the POSCAR_moved directory does not contain previous data
        if os.path.exists(moved_dir):
            shutil.rmtree(moved_dir)

        # Move monolayer structure to the specified interlayer position
        move_poscar(monolayer_file, position)

        # Retrieve the newly moved monolayer structure file
        moved_monolayer_file = os.path.join(moved_dir, os.listdir(moved_dir)[0])

        # Construct the bilayer structure and write the POSCAR files saved in the nested directory
        Bilayer(moved_monolayer_file, monolayer_file, overwrite=False, skip_xy_rev=True).WritePOSCAR()
        # Optionally, all the POSCAR files are saved in the same folder f"{savepath}/{formula}_{genmode}_{cord*}-{mismatch}-{natom1}-POSCAR"
        # Bilayer(moved_monolayer_file, monolayer_file, overwrite=False, skip_xy_rev=True, savenamemode=2).WritePOSCAR()

shutil.rmtree(moved_dir)

# After using Bilayer with savenamemode=2, rename the POSCAR filenames with ".vasp"
# for file in os.listdir("BiPOSCAR_dir"):
#     shutil.move(os.path.join("BiPOSCAR_dir", file), os.path.join("BiPOSCAR_dir", file.replace("-POSCAR", ".vasp")))
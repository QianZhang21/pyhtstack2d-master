from pymatgen.core import Structure
import numpy as np
import os

# Read structure file paths and magnetic moments from input file
with open("mid.txt", "r") as f:
    lines = f.readlines()

mid_flat = []

# Loop through each line in mid.txt
for li in lines:
    filepath = li.strip().split(" ")[0]  # Extract path
    mag = li.strip().split(" ")[1]       # Extract magnetic moment

    # Filter out structures with magnetic moment > 0.1
    if abs(float(mag)) > 0.1:
        continue

    # Construct full path to the structure file (CONTCAR)
    stfile = os.path.join(filepath.strip(), "AA", "cord1", "CONTCAR")

    # Load the structure using pymatgen
    structure = Structure.from_file(stfile)
    labels = structure.labels
    positions = structure.cart_coords

    B_z_positions = []  # z-coordinates of B atoms
    N_z_positions = []  # z-coordinates of N atoms
    num_BN = 0          # Count of B atoms (equal to N atoms assumed)

    # Separate B and N atom z-coordinates
    for i, label in enumerate(labels):
        if label == "B":
            num_BN += 1
            B_z_positions.append(positions[i][2])
        elif label == "N":
            N_z_positions.append(positions[i][2])

    # Sort z-coordinates for B and N atoms
    B_z_positions.sort(key=lambda x: x)
    N_z_positions.sort(key=lambda x: x)

    # Divide into top and bottom BN layer atoms (assuming bilayer)
    BN_up = B_z_positions[:int(num_BN/2)] + N_z_positions[:int(num_BN/2)]
    BN_down = B_z_positions[int(num_BN/2):] + N_z_positions[int(num_BN/2):]
    BN_up.sort(key=lambda x: x)
    BN_down.sort(key=lambda x: x)

    # Check if both BN layers are "flat" (i.e., minimal z spread)
    if (np.max(BN_up) - np.min(BN_up)) < 0.1 and (np.max(BN_down) - np.min(BN_down)) < 0.1:
        mid_flat.append(filepath.strip())
        # Uncomment for debugging:
        # print(filepath.strip(), f"{np.max(BN_up)-np.min(BN_up)}, {np.max(BN_down)-np.min(BN_down)}")
    else:
        print(filepath.strip(), f"{np.max(BN_up)-np.min(BN_up)}, {np.max(BN_down)-np.min(BN_down)}")

# Write paths of flat bilayer BN structures to a file
with open("mid_flat.txt", "w", newline="\n") as f:
    for item in mid_flat:
        f.write("%s\n" % item)



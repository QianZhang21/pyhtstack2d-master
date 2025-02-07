import os

import numpy as np


def rotationmatrix(angle, axis="z"):
    """
    Generate a rotation matrix for a given angle (degree) and axis (x, y, z)

    angle: float
        The angle of rotation in degree
    axis: str
        The axis of rotation
    """
    angle = np.deg2rad(angle)
    if axis.lower() == "z":
        return np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0],
                         [0, 0, 1]])
    elif axis.lower() == "y":
        return np.array([[np.cos(angle), 0, np.sin(angle)],
                         [0, 1, 0],
                         [-np.sin(angle), 0, np.cos(angle)]])
    elif axis.lower() == "x":
        return np.array([[1, 0, 0],
                         [0, np.cos(angle), -np.sin(angle)],
                         [0, np.sin(angle), np.cos(angle)]])
    else:
        raise ValueError("Invalid axis")


def read_poscar(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    title = lines[0].strip()
    scaling_factor = float(lines[1].strip())
    lattice_vectors = np.array([list(map(float, line.strip().split())) for line in lines[2:5]])
    atom_types = lines[5].strip().split()
    atom_counts = list(map(int, lines[6].strip().split()))
    coord_type = lines[7].strip().lower()
    atom_positions = np.array([list(map(float, line.strip().split()[:3])) for line in lines[8:8 + sum(atom_counts)]])
    return title, scaling_factor, lattice_vectors, atom_types, atom_counts, coord_type, atom_positions


def write_poscar(file_path, title, scaling_factor, lattice_vectors, atom_types, atom_counts, coord_type, atom_positions):
    with open(file_path, 'w') as file:
        file.write(f"{title}\n")
        file.write(f"{scaling_factor}\n")
        for vector in lattice_vectors:
            file.write(" ".join(f"{x:.16f}" for x in vector) + "\n")
        file.write(" ".join(atom_types) + "\n")
        file.write(" ".join(map(str, atom_counts)) + "\n")
        file.write(f"{coord_type.capitalize()}\n")
        for position in atom_positions:
            file.write(" ".join(f"{x:.16f}" for x in position) + "\n")


def rotate_poscar(file_path, angle, output_path=None):
    if output_path is None:
        output_count = 0
        if not os.path.exists("POSCAR_rotated"):
            os.makedirs("POSCAR_rotated")
        output_path = os.path.join("POSCAR_rotated", os.path.basename(file_path))
        while os.path.exists(output_path):
            output_count += 1
            output_path = os.path.join("POSCAR_rotated", f"{os.path.basename(file_path)}_{output_count}")
    title, scaling_factor, lattice_vectors, atom_types, atom_counts, coord_type, atom_positions = read_poscar(file_path)
    # Rotate only the atomic positions while preserving the lattice vectors
    # Generate the rotation matrix for a rotation about the z-axis
    rotation_matrix = rotationmatrix(angle, axis='z')
    # Convert atomic positions to Cartesian coordinates
    # Convert to Cartesian if coordinates are in Direct format
    if coord_type.lower().startswith('d'):
        cartesian_positions = np.dot(atom_positions, lattice_vectors)
        rotated_cartesian_positions = np.dot(cartesian_positions, rotation_matrix.T)
        rotated_positions = np.dot(rotated_cartesian_positions, np.linalg.inv(lattice_vectors))
    else:
        rotated_positions = np.dot(atom_positions, rotation_matrix.T)
    # Write the new POSCAR with rotated atomic positions
    write_poscar(output_path, title, scaling_factor, lattice_vectors, atom_types, atom_counts, coord_type, rotated_positions)
    print(f"Rotated POSCAR with preserved structure saved to {output_path}")


def move_poscar(file_path, translation_vector, output_path=None):
    if output_path is None:
        output_count = 0
        if not os.path.exists("POSCAR_moved"):
            os.makedirs("POSCAR_moved")
        output_path = os.path.join("POSCAR_moved", os.path.basename(file_path))
        while os.path.exists(output_path):
            output_count += 1
            output_path = os.path.join("POSCAR_moved", f"{os.path.basename(file_path)}_{output_count}")
    title, scaling_factor, lattice_vectors, atom_types, atom_counts, coord_type, atom_positions = read_poscar(file_path)
    # Move the atomic positions by the translation vector
    if np.all(np.abs(translation_vector) <= 1) and coord_type.lower().startswith('c'):
        fractions_positions = np.dot(atom_positions, np.linalg.inv(lattice_vectors))
        moved_positions = fractions_positions + translation_vector
        # Wrap positions to ensure they are within [0, 1)
        wrapped_positions = moved_positions % 1.0
        wrapped_positions = np.dot(wrapped_positions, lattice_vectors)
    else:
        moved_positions = atom_positions + translation_vector
        # Wrap positions to ensure they are within [0, 1)
        if coord_type.lower().startswith('d'):
            wrapped_positions = moved_positions % 1.0
        else:
            wrapped_positions = np.dot(moved_positions, np.linalg.inv(lattice_vectors))
            wrapped_positions = wrapped_positions % 1.0
            wrapped_positions = np.dot(wrapped_positions, lattice_vectors)
    write_poscar(output_path, title, scaling_factor, lattice_vectors, atom_types, atom_counts, coord_type, wrapped_positions)
    print(f"Modified POSCAR with positions moved and wrapped within the unit cell saved to {output_path}")
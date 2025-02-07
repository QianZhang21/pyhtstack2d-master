import numpy as np
from pymatgen.io.vasp import Poscar
from pymatgen.transformations.advanced_transformations import SupercellTransformation
from pymatgen.core import Lattice
import os


def transformlattice(filename="POSCAR", transmat=None, overwrite=False):
    """
    Transform the lattice of a POSCAR file.

    filename: str
        The name of the POSCAR file.
    transmat: numpy.ndarray
        The transformation matrix. If None, the default matrix is used.
        For example, transmat = np.array([[-1, 1, 0], [-1, -1, 0], [0, 0, 1]]) will transform the hexagonal lattice to a orthorhombic lattice.
    overwrite: bool
        Whether to overwrite the original POSCAR file.
    """
    assert transmat is None or (isinstance(transmat, np.ndarray) and transmat.shape == (3, 3)), \
        "transmat must be a 3x3 numpy array or None"

    if transmat is None:
        transmat = np.array([[-1, 1, 0], [-1, -1, 0], [0, 0, 1]])

    poscar = Poscar.from_file(filename)
    transformation = SupercellTransformation(scaling_matrix=transmat.tolist())

    transformst = transformation.apply_transformation(poscar.structure)
    transformla = transformst.lattice.matrix

    cosab = np.dot(transformla[0], transformla[1]) / (np.linalg.norm(transformla[0]) * np.linalg.norm(transformla[1]))
    cosac = np.dot(transformla[0], transformla[2]) / (np.linalg.norm(transformla[0]) * np.linalg.norm(transformla[2]))
    laa = [np.linalg.norm(transformla[0]), 0, 0]
    lab = [np.linalg.norm(transformla[1]) * cosab, np.linalg.norm(transformla[1]) * np.sqrt(1 - cosab ** 2), 0]
    lac = [np.linalg.norm(transformla[2]) * cosac, 0, np.linalg.norm(transformla[2]) * np.sqrt(1 - cosac ** 2)]
    transformst.lattice = Lattice([laa, lab, lac])

    transposcar = Poscar(transformst)
    if overwrite:
        savefilename = filename
    else:
        count = 1
        savefilename = filename + "_transformed"
        while os.path.exists(savefilename):
            savefilename = filename + "_transformed_" + str(count)
            count += 1
    transposcar.write_file(savefilename)
    print("The transformed POSCAR file is saved as", savefilename)



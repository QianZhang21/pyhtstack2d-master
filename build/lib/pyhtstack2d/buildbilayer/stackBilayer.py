import os
import warnings
from collections import Counter

import numpy as np

rotation_matrix_dict = {0: np.array([[1., 0., 0.],
                                     [0., 1., 0.],
                                     [0., 0., 1.]]),
                        60: np.array([[0.5, -0.8660254, 0.],
                                      [0.8660254, 0.5, 0.],
                                      [0., 0., 1.]]),
                        90: np.array([[0., -1., 0.],
                                      [1., 0., 0.],
                                      [0., 0., 1.]]),
                        120: np.array([[-0.5, -0.8660254, 0.],
                                       [0.8660254, -0.5, 0.],
                                       [0., 0., 1.]]),
                        180: np.array([[-1., 0., 0.],
                                       [0., -1., 0.],
                                       [0., 0., 1.]])}


def match_positions(pos1, pos2, lattice_vectors, angles=None):
    # Remove duplicate structures
    # Rotate pos2 using the provided rotation function
    if angles is None:
        angles = [0, 60, 90, 120, 180]
    for angle in angles:
        # Rotate pos2 using the provided rotation function
        # Generate the rotation matrix for a rotation about the z-axis
        rotation_matrix = rotation_matrix_dict[angle]
        # Convert atomic positions to Cartesian coordinates
        cartesian_positions = np.dot(pos2, lattice_vectors)
        # Rotate the atomic positions in Cartesian space
        rotated_cartesian_positions = np.dot(cartesian_positions, rotation_matrix.T)
        # Convert back to fractional coordinates and normalize
        rotated_pos2 = np.round(np.dot(rotated_cartesian_positions, np.linalg.inv(lattice_vectors)) % 1.0, 6)
        if np.allclose(rotated_pos2, pos1, atol=1e-3):
            return True
    return False


class Monolayer:
    def __init__(self, stfile):
        """
        Get the 2D structure  from a POSCAR file.

        stfile: str
            Path to the POSCAR file.
            Default Z-axis is vacuum layer (out-of-plane) direction.
        """
        self.StFile = stfile
        self.uid = os.path.basename(self.StFile).split('-')[0]
        self.lattice = None
        self.a = None
        self.b = None
        self.c = None
        self.cos_angle1 = None
        self.cos_angle2 = None
        self.cos_angle3 = None
        self.elements = None
        self.elemnum = None
        self.posType = None
        self.positions = None
        self.positions_direct = None
        self.tmdtype = None
        self.getInfo()

    def getInfo(self):
        """
        Get the information from the structure file, ensuring it exists and parsing it.
        """
        if not os.path.exists(self.StFile):
            raise FileNotFoundError("The POSCAR file does not exist.")

        with open(self.StFile, "r") as f:
            lines = f.readlines()

        self.lattice = np.array([list(map(float, lines[i].split())) for i in range(2, 5)])
        self.elements = list(map(str, lines[5].split()))
        self.elemnum = list(map(int, lines[6].split()))
        self.posType = lines[7].strip().lower()
        positions_ = np.array([list(map(float, line.split()[:3])) for line in lines[8:8 + sum(self.elemnum)]])
        assert self.lattice.shape == (3, 3), "Lattice vectors are not 3x3."
        assert len(self.elements) == len(self.elemnum), "Number of elements and element numbers do not match."

        self.a = np.linalg.norm(self.lattice[0])
        self.b = np.linalg.norm(self.lattice[1])
        # if np.linalg.norm(self.lattice[0]) >= np.linalg.norm(self.lattice[1]):
        #     self.a = np.linalg.norm(self.lattice[0])
        #     self.b = np.linalg.norm(self.lattice[1])
        # else:
        #     self.a = np.linalg.norm(self.lattice[1])
        #     self.b = np.linalg.norm(self.lattice[0])
        self.c = np.linalg.norm(self.lattice[2])
        self.cos_angle1 = np.dot(self.lattice[0], self.lattice[1]) / (self.a * self.b)
        self.cos_angle2 = np.dot(self.lattice[0], self.lattice[2]) / (self.a * self.c)
        self.cos_angle3 = np.dot(self.lattice[1], self.lattice[2]) / (self.b * self.c)

        if self.a > self.c or self.b > self.c:
            raise ValueError("Please check the vacuum layer direction: z-axis is not the longest lattice vector.")

        if self.posType.startswith('d'):
            self.positions_direct = positions_.copy()
            self.positions = np.dot(positions_, self.lattice)
        elif self.posType.startswith('c'):
            self.positions = positions_.copy()
            self.positions_direct = np.dot(positions_, np.linalg.inv(self.lattice))
        else:
            raise ValueError("Please check the position type: it is not cartesian or direct.")

        # Sorts positions and corresponding elements by the z-coordinate.
        # Get indices that would sort the positions array by the z-coordinate (3rd column)
        indices = np.argsort(self.positions[:, 2])
        self.positions = self.positions[indices]
        self.positions_direct = self.positions_direct[indices]
        sorted_elements = []
        for elem, num in zip(self.elements, self.elemnum):
            sorted_elements.extend([elem] * num)
        self.elements = [sorted_elements[i] for i in indices]

    def get_intra_distances(self):
        """
        Calculates the differences in z-coordinates between consecutive atoms.
        Returns a numpy array of these differences.
        """
        z_coords = self.positions[:, 2]  # Extract just the z-coordinates
        return np.diff(z_coords)

    def checkLatticeType(self):
        """
        Checks the type of lattice based on the lattice vectors and angles for the 2D monolayer.
        The accuracy of crystal system detection can be improved by changing the atol.
        """
        atol = 0.1
        if np.abs(self.a - self.b) < atol:
            if np.isclose(self.cos_angle1, 0, atol=atol):
                return "Square"
            elif np.isclose(self.cos_angle1, 0.5, atol=atol) or np.isclose(self.cos_angle1, -0.5, atol=atol):
                return "Hexagonal"
            else:
                return None
        else:
            if np.isclose(self.cos_angle1, 0, atol=atol):
                return "Rectangular"
            else:
                return "Oblique"

    def checkTMDH(self, natoms=3):
        """
        When the material is a TMD/MH/TMO, it is checked if it is a 1T or 2H type, and the lattice constant is returned.
        Checking for a hexagonal lattice by comparing angles and norms of lattice vectors.
        """
        if sum(self.elemnum) != natoms:
            raise ValueError("Please check if material is a TMD/MH/TMO: Number of elements is not 3.")

        if not (np.isclose(self.cos_angle1, 0.5, atol=0.01) or np.isclose(self.cos_angle1, -0.5, atol=0.01)):
            raise ValueError("Please check if lattice is hexagonal: Lattice angles are not 60 or 120 degrees.")

        if not (np.isclose(self.cos_angle2, 0, atol=0.01) and np.isclose(self.cos_angle3, 0, atol=0.01)):
            raise ValueError("Please check if lattice is hexagonal: Lattice angles with C are not 90 degrees.")

        if abs(self.a - self.b) > 0.01:
            raise ValueError("Please check if lattice is hexagonal: Lattice vectors a and b are not equal.")

        if natoms == 3:
            if np.linalg.norm(self.positions[0, :2] - self.positions[2, :2]) < 0.1:
                self.tmdtype = '2H'
            else:
                self.tmdtype = '1T'
        elif natoms == 4:
            if np.linalg.norm(self.positions[0, :2] - self.positions[3, :2]) < 0.1 and \
                    np.linalg.norm(self.positions[1, :2] - self.positions[2, :2]) < 0.1:
                self.tmdtype = '2H'
            else:
                self.tmdtype = '1T'
        else:
            raise ValueError("Please check the POSCAR: Number of elements is not 3 or 4.")
        return self.a


class Bilayer:
    def __init__(self, st1, st2=None, la1=None, la2=None, position1=None, position2=None, d_inter=3.3, lv=30.0,
                 savepath=None, savenamemode=1, hetero=True, mismatch_threshold=5.0, stackmode="AA", overwrite=True,
                 skip_xy_rev=False):
        """
        Constructing 2D bilayer structures.

        st1: list or str
            The atoms in the first layer or the POSCAR file of the first layer.
            If st1 is a list, note that the sequence of atoms in the list is ordered by height in Z direction.
            If st1 is a str, it is the path to the POSCAR file of the first layer, and a1, d_intra1 will be ignored.
        st2: list or str
            The atoms in the second layer or the POSCAR file of the second layer.
            If not provided, homogeneous bilayer will be constructed.
        la1: array
            The in-plane lattice constant for the first layer, which has the up_shape (3, 3).
        la2: array
            The in-plane lattice constant for the second layer.
            If not provided, it will be set to a1.
        position1: array
            The positions of atoms in the first layer, which has the up_shape (n, 3).
            Please provide the positions in the direct coordinates.
        position2: array
            The positions of atoms in the second layer.
            Please provide the positions in the direct coordinates.
        d_inter: float
            The distance between two atoms in different layers.
        lv: float
            The length of the vacuum layer.
        savepath: str
            The path to save the POSCAR files.
            If not provided, will be saved in the BiPOSCAR_dir directory.
        savenamemode: int
            The mode to name the saved POSCAR files.
            1: The POSCAR files are saved in the folder f"{savepath}/{formula}_{mismatch}_{natom1}/{genmode}/{cord*}/POSCAR".
            2: All the POSCAR files are saved in the same folder f"{savepath}/{formula}_{genmode}_{cord*}-{mismatch}-{natom1}-POSCAR".
            (mismatch is the mismatch between two layers, natom1 is the number of atoms in the first layer, cord* is the different stacking sequences.)
        hetero: bool
            If True, regardless of whether the elements of st1 and st2 are consistent, they are stacked as heterogeneous bilayers, mainly affecting the atomic arrangement in the z-axis direction.
            If False, they stack as homogeneous bilayers when the elements are consistent, e.g. when hetero=False, the 1T and 2H MoS2 bilayers fail to stack.
        mismatch_threshold: float
            Maximum tolerance of mismatch_threshold lattice constant between two layers.
        genmode: str
             The parameter in this class is only used as a folder name for saving.
        stackmode: str
            The stacking mode name of the bilayer
        overwrite: bool
            If True, the existing POSCAR files will be overwritten.
        skip_xy_rev: bool
            When skip_xy_rev is true, it does not take into account whether the XY position can be rotated. For example, for a 1T TMD monolayer, a 60° rotation will not be considered.
        """
        if savepath is None:
            self.savepath = 'BiPOSCAR_dir'
        else:
            self.savepath = savepath
        if not os.path.exists(self.savepath):
            os.makedirs(self.savepath)
        if savenamemode not in [1, 2]:
            warnings.warn("The savenamemode is not valid. It will be set to 1.")
            savenamemode = 1
        self.savenamemode = savenamemode
        self.mismatch_threshold = mismatch_threshold
        self.mismatch = 0.0
        self.skip_xy_rev = skip_xy_rev
        self.overwrite = overwrite

        if isinstance(st1, str) or isinstance(st1, Monolayer):
            if isinstance(st1, str):
                mono1 = Monolayer(st1)
            else:
                mono1 = st1
            Latype1 = mono1.checkLatticeType()
            a1 = mono1.a
            b1 = mono1.b
            angle1 = mono1.cos_angle1
            elements1 = mono1.elements
            pos1_xy = mono1.positions_direct[:, :2]
            pos1_z = mono1.positions[:, 2]
            lattice_w = mono1.lattice
            lattice_w[2, 2] = lv
            indices1 = None
        elif isinstance(st1, list):
            assert la1 is not None and position1 is not None, "Please provide the in-plane lattice constant and the " \
                                                              "positions for the first layer."
            assert len(st1) == len(position1), "The number of atoms in the first layer and the number of positions " \
                                               "provided are not the same."
            assert np.array(la1).shape == (3, 3), "Please provide the lattice constant of the first layer in the " \
                                                  "up_shape of (3, 3)."
            Latype1, a1, b1, angle1 = self.checkLatticeType(np.array(la1))
            assert np.max(np.array(position1)) < 1, "Please provide the positions in the direct coordinates."
            elements1, pos1_xy, pos1_z, indices1 = self.sorted_position(st1, np.array(position1), np.array(la1))
            lattice_w = np.r_[np.array(la1)[:2], [[0, 0, lv]]]
        else:
            raise ValueError("Please provide the atoms in the first layer or the POSCAR file of the first layer.")

        w_mode = -1
        if (isinstance(st2, str) or isinstance(st2, Monolayer)) and st1 != st2:
            if isinstance(st2, str):
                mono2 = Monolayer(st2)
            else:
                mono2 = st2
            elements2 = mono2.elements
            if (elements1 == elements2 or elements1 == elements2[::-1]) and not hetero:
                w_mode = 1
                elements2 = elements1
                Latype2, a2, b2, angle2, pos2_xy, pos2_z = Latype1, a1, b1, angle1, pos1_xy, pos1_z
            else:
                Latype2 = mono2.checkLatticeType()
                a2 = mono2.a
                b2 = mono2.b
                angle2 = mono2.cos_angle1
                pos2_xy = mono2.positions_direct[:, :2]
                pos2_z = mono2.positions[:, 2]
        elif isinstance(st2, list):
            if (st2 == elements1 or st2 == elements1[::-1]) and not hetero:
                w_mode = 1
                elements2 = elements1
                Latype2, a2, b2, angle2, pos2_xy, pos2_z = Latype1, a1, b1, angle1, pos1_xy, pos1_z
            else:
                if la2 is not None:
                    assert np.array(la2).shape == (3, 3), "Please provide the lattice constant of the second layer in " \
                                                          "the up_shape of (3, 3)."
                    assert position2 is not None, "Please provide the positions for the second layer."
                    assert len(st2) == len(position2), "The number of atoms in the second layer and the number of " \
                                                       "positions provided are not the same."
                    Latype2, a2, b2, angle2 = self.checkLatticeType(np.array(la2))
                    assert np.max(np.array(position2)) < 1, "Please provide the positions in the direct coordinates."
                    elements2, pos2_xy, pos2_z, indices2 = self.sorted_position(st2, np.array(position2), np.array(la2))
                else:
                    if indices1 is not None:
                        elements2 = [st2[i] for i in indices1]
                    else:
                        elements2 = st2
                    Latype2, a2, b2, angle2, pos2_xy, pos2_z = Latype1, a1, b1, angle1, pos1_xy, pos1_z
        else:
            w_mode = 1
            Latype2, a2, b2, angle2, elements2, pos2_xy, pos2_z = Latype1, a1, b1, angle1, elements1, pos1_xy, pos1_z

        assert Latype1 == Latype2, "The lattice type of the two layers are different."

        self.Latype = Latype1
        self.lv = lv
        self.stackmode = stackmode

        if not hetero:
            self.elem1_, self.elem2_ = self.reduce_element_pos(elements1, elements2, pos1_z, pos2_z)
        else:
            self.elem1_, self.elem2_ = elements1, elements2

        self.natom1 = len(elements1)

        if w_mode == 1:
            self.hetero = False
            self.lattice_w = lattice_w
            self.elements = elements1 + elements1
            if self.elem1_ == self.elem1_[::-1]:
                self.seq_cord = ["cord1"]
            else:
                self.seq_cord = ["cord1", "cord2", "cord3"]
        else:
            self.hetero = True
            self.lattice_w = self.check_lattice_mismatch(a1, b1, angle1, a2, b2, angle2)
            self.elements = elements1 + elements2
            if self.elem1_ == self.elem1_[::-1] and self.elem2_ == self.elem2_[::-1]:
                self.seq_cord = ["cord1"]
            elif self.elem1_ == self.elem1_[::-1] and self.elem2_ != self.elem2_[::-1]:
                self.seq_cord = ["cord1", "cord3"]
            elif self.elem1_ != self.elem1_[::-1] and self.elem2_ == self.elem2_[::-1]:
                self.seq_cord = ["cord1", "cord2"]
            else:
                self.seq_cord = ["cord1", "cord2", "cord3", "cord4"]
        count = Counter(self.elements)
        self.formula_w = ''.join(f"{element}{count[element] if count[element] > 1 else ''}" for element in count)
        if self.savenamemode == 1:
            self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + f"{self.natom1}"

        self.d_inter = d_inter
        self.pos1_xy = pos1_xy
        self.pos2_xy = pos2_xy
        self.pos1_z = pos1_z
        self.pos2_z = pos2_z
        self.ini_position()

    def checkLatticeType(self, la):
        """
        Checks the type of lattice based on the lattice vectors and angles for the 2D monolayer.
        The accuracy of crystal system detection can be improved by changing the atol.
        """
        atol = 0.1
        if np.linalg.norm(la[0]) >= np.linalg.norm(la[1]):
            a = np.linalg.norm(la[0])
            b = np.linalg.norm(la[1])
        else:
            a = np.linalg.norm(la[1])
            b = np.linalg.norm(la[0])
        cos_angle1 = np.dot(la[0], la[1]) / (a * b)
        if np.abs(a - b) < atol:
            if np.isclose(cos_angle1, 0, atol=atol):
                return "Square", a, b, cos_angle1
            elif np.isclose(cos_angle1, 0.5, atol=atol) or np.isclose(cos_angle1, -0.5, atol=atol):
                return "Hexagonal", a, b, cos_angle1
            else:
                return None, a, b, cos_angle1
        else:
            if np.isclose(cos_angle1, 0, atol=atol):
                return "Rectangular", a, b, cos_angle1
            else:
                return "Oblique", a, b, cos_angle1

    def check_lattice_mismatch(self, la1, lb1, angle1, la2, lb2, angle2):
        assert self.Latype, "The lattice type is not supported."
        mislattice_a = 2 * abs(la1 - la2) / (la1 + la2) * 100
        mislattice_b = 2 * abs(lb1 - lb2) / (lb1 + lb2) * 100
        self.mismatch = max(mislattice_a, mislattice_b)
        la_mean = (la1 + la2) / 2
        lb_mean = (lb1 + lb2) / 2
        if mislattice_a > self.mismatch_threshold:
            warnings.warn(
                f"The x-axis lattice mismatch_threshold {mislattice_a:.2f}% is larger than {self.mismatch_threshold}%.")
        if self.Latype == "Square":
            return np.array([[la_mean, 0, 0], [0, la_mean, 0], [0, 0, self.lv]])
        elif self.Latype == "Hexagonal":
            return np.array(
                [[la_mean, 0, 0], [la_mean * angle1, la_mean * np.sqrt(1 - angle1 ** 2), 0], [0, 0, self.lv]])
        elif self.Latype == "Rectangular":
            if mislattice_b > self.mismatch_threshold:
                warnings.warn(
                    f"The y-axis lattice mismatch_threshold {mislattice_b:.2f}% is larger than {self.mismatch_threshold}%.")
            return np.array([[la_mean, 0, 0], [0, lb_mean, 0], [0, 0, self.lv]])
        elif self.Latype == "Oblique" and np.isclose(angle1, angle2, atol=0.1):
            if mislattice_b > self.mismatch_threshold:
                warnings.warn(
                    f"The y-axis lattice mismatch_threshold {mislattice_b:.2f}% is larger than {self.mismatch_threshold}%.")
            return np.array(
                [[la_mean, 0, 0], [lb_mean * angle1, lb_mean * np.sqrt(1 - angle1 ** 2), 0], [0, 0, self.lv]])
        else:
            raise ValueError("Please check the lattice vectors and angles for the two layers.")

    def reduce_element_pos(self, elements1, elements2, pos1_z, pos2_z):
        """
        Reduces the elements and positions to the elements that are not in the same plane.
        After reduction, some stacking orders may be ignored.
        """
        elem1 = []
        elem2 = []
        z_coords_diff1 = np.diff(pos1_z)
        z_coords_diff2 = np.diff(pos2_z)
        for i, z_diff in enumerate(z_coords_diff1):
            if z_diff > 0.01:
                if len(elem1) == 0:
                    elem1.append(elements1[i])
                    elem1.append(elements1[i + 1])
                else:
                    elem1.append(elements1[i + 1])

        for i, z_diff in enumerate(z_coords_diff2):
            if z_diff > 0.01:
                if len(elem2) == 0:
                    elem2.append(elements2[i])
                    elem2.append(elements2[i + 1])
                else:
                    elem2.append(elements2[i + 1])
        return elem1, elem2

    def sorted_position(self, elem, pos, la):
        """
        Orders the positions based on the height in Z direction.
        """
        pos_c = np.dot(pos, la)
        indices = np.argsort(pos[:, 2])
        sorted_pos = pos[indices]
        sorted_pos_c = pos_c[indices]
        sorted_elem = [elem[i] for i in indices]
        return sorted_elem, sorted_pos, sorted_pos_c[:, 2], indices

    def ini_position(self):
        min_pos1 = np.min(self.pos1_z)
        max_pos1 = np.max(self.pos1_z)
        min_pos2 = np.min(self.pos2_z)
        max_pos2 = np.max(self.pos2_z)
        dis1 = max_pos1 - min_pos1
        dis2 = max_pos2 - min_pos2

        if dis1 * 2 < dis2:
            self.pos1_z = self.pos1_z - min_pos1 + self.d_inter + self.lv / 2
            self.pos2_z = self.pos2_z - min_pos2 + self.lv / 2 - (max_pos2 - min_pos2)
        elif dis2 * 2 < dis1:
            self.pos1_z = self.pos1_z - min_pos1 + self.lv / 2
            self.pos2_z = self.pos2_z - min_pos2 + self.lv / 2 - (max_pos2 - min_pos2) - self.d_inter
        else:
            self.pos1_z = self.pos1_z - min_pos1 + self.d_inter / 2 + self.lv / 2
            self.pos2_z = self.pos2_z - min_pos2 + self.lv / 2 - (max_pos2 - min_pos2) - self.d_inter / 2
        self.pos1_z = self.pos1_z / self.lv
        self.pos2_z = self.pos2_z / self.lv

        pos1_shift = np.column_stack((self.pos1_xy[:, 0], self.pos1_xy[:, 1], self.pos1_z))
        pos2_shift = np.column_stack((self.pos2_xy[:, 0], self.pos2_xy[:, 1], self.pos2_z))

        pos1_z_copy = self.pos1_z.copy()[::-1]
        pos2_z_copy = self.pos2_z.copy()[::-1]
        d_intra1 = np.diff(self.pos1_z)
        d_intra2 = np.diff(self.pos2_z)

        for i, d1 in enumerate(d_intra1):
            pos1_z_copy[i + 1] = pos1_z_copy[i] - d1
        for i, d2 in enumerate(d_intra2):
            pos2_z_copy[i + 1] = pos2_z_copy[i] - d2

        rev_pos1_shift = np.column_stack((self.pos1_xy[:, 0], self.pos1_xy[:, 1], pos1_z_copy))
        rev_pos2_shift = np.column_stack((self.pos2_xy[:, 0], self.pos2_xy[:, 1], pos2_z_copy))

        self.pos_w = {self.stackmode: {"cord1": np.concatenate((pos1_shift, pos2_shift), axis=0),
                                       "cord2": np.concatenate((rev_pos1_shift, pos2_shift), axis=0),
                                       "cord3": np.concatenate((pos1_shift, rev_pos2_shift), axis=0),
                                       "cord4": np.concatenate((rev_pos1_shift, rev_pos2_shift), axis=0)}}

        if not self.skip_xy_rev:
            pos1_xy_rev, pos2_xy_rev = False, False
            # pos1_xy_copy, pos2_xy_copy = None, None
            if self.elem1_ == self.elem1_[::-1]:
                pos1_xy_copy = self.pos1_xy.copy()[::-1]
                for pos1_xy_i in range(int(len(self.pos1_xy) / 2)):
                    if np.linalg.norm(pos1_xy_copy[pos1_xy_i] - self.pos1_xy[pos1_xy_i]) > 0.01:
                        pos1_xy_rev = True
                        break
                if pos1_xy_rev:
                    pos1_shift = np.column_stack((pos1_xy_copy[:, 0], pos1_xy_copy[:, 1], self.pos1_z))
                    # pos2_shift = np.column_stack((self.pos2_xy[:, 0], self.pos2_xy[:, 1], self.pos2_z))
                    rev_pos1_shift = np.column_stack((pos1_xy_copy[:, 0], pos1_xy_copy[:, 1], pos1_z_copy))
                    # rev_pos2_shift = np.column_stack((self.pos2_xy[:, 0], self.pos2_xy[:, 1], pos2_z_copy))
                    self.pos_w[self.stackmode + str(2)] = {"cord1": np.concatenate((pos1_shift, pos2_shift), axis=0),
                                                           "cord2": np.concatenate((rev_pos1_shift, pos2_shift),
                                                                                   axis=0),
                                                           "cord3": np.concatenate((pos1_shift, rev_pos2_shift),
                                                                                   axis=0),
                                                           "cord4": np.concatenate((rev_pos1_shift, rev_pos2_shift),
                                                                                   axis=0)}
            stackmode3 = True
            if np.linalg.norm(self.pos2_xy[0] - self.pos1_xy[0]) < 0.01 or np.linalg.norm(
                    self.pos2_xy[0] - self.pos1_xy[-1]) < 0.01:
                if np.linalg.norm(self.pos2_xy[-1] - self.pos1_xy[0]) < 0.01 or np.linalg.norm(
                        self.pos2_xy[-1] - self.pos1_xy[-1]) < 0.01:
                    stackmode3 = False

            if self.elem2_ == self.elem2_[::-1] and stackmode3:
                pos2_xy_copy = self.pos2_xy.copy()[::-1]
                for pos2_xy_i in range(int(len(self.pos2_xy) / 2)):
                    if np.linalg.norm(pos2_xy_copy[pos2_xy_i] - self.pos2_xy[pos2_xy_i]) > 0.01:
                        pos2_xy_rev = True
                        break
                if self.hetero and pos2_xy_rev:
                    pos1_shift = np.column_stack((self.pos1_xy[:, 0], self.pos1_xy[:, 1], self.pos1_z))
                    pos2_shift = np.column_stack((pos2_xy_copy[:, 0], pos2_xy_copy[:, 1], self.pos2_z))
                    rev_pos1_shift = np.column_stack((self.pos1_xy[:, 0], self.pos1_xy[:, 1], pos1_z_copy))
                    rev_pos2_shift = np.column_stack((pos2_xy_copy[:, 0], pos2_xy_copy[:, 1], pos2_z_copy))
                    self.pos_w[self.stackmode + str(3)] = {"cord1": np.concatenate((pos1_shift, pos2_shift), axis=0),
                                                           "cord2": np.concatenate((rev_pos1_shift, pos2_shift),
                                                                                   axis=0),
                                                           "cord3": np.concatenate((pos1_shift, rev_pos2_shift),
                                                                                   axis=0),
                                                           "cord4": np.concatenate((rev_pos1_shift, rev_pos2_shift),
                                                                                   axis=0)}
                    if self.stackmode + str(2) in self.pos_w:
                        if match_positions(self.pos_w[self.stackmode + str(2)]["cord1"],
                                           self.pos_w[self.stackmode + str(3)]["cord1"], self.lattice_w):
                            self.pos_w.pop(self.stackmode + str(3))

            # if self.hetero and pos1_xy_rev and pos2_xy_rev:
            #     pos1_shift = np.column_stack((pos1_xy_copy[:, 0], pos1_xy_copy[:, 1], self.pos1_z))
            #     pos2_shift = np.column_stack((pos2_xy_copy[:, 0], pos2_xy_copy[:, 1], self.pos2_z))
            #     rev_pos1_shift = np.column_stack((pos1_xy_copy[:, 0], pos1_xy_copy[:, 1], pos1_z_copy))
            #     rev_pos2_shift = np.column_stack((pos2_xy_copy[:, 0], pos2_xy_copy[:, 1], pos2_z_copy))
            #     self.pos_w[self.stackmode + str(4)] = {"cord1": np.concatenate((pos1_shift, pos2_shift), axis=0),
            #                                            "cord2": np.concatenate((rev_pos1_shift, pos2_shift), axis=0),
            #                                            "cord3": np.concatenate((pos1_shift, rev_pos2_shift), axis=0),
            #                                            "cord4": np.concatenate((rev_pos1_shift, rev_pos2_shift),
            #                                                                    axis=0)}

    def WritePOSCAR(self):
        """
        Writes the POSCAR file for the bilayer.
        """
        for stackmode_ in self.pos_w.keys():
            for cord_ in self.seq_cord:
                element_pos_dict = {}
                for i, e in enumerate(self.elements):
                    if e not in element_pos_dict:
                        element_pos_dict[e] = {}
                        element_pos_dict[e]['count'] = 0
                        element_pos_dict[e]['pos'] = []
                    element_pos_dict[e]['count'] += 1
                    element_pos_dict[e]['pos'].append(self.pos_w[stackmode_][cord_][i])

                if self.savenamemode == 1:
                    stackmode_add_count = stackmode_
                    if not os.path.exists(os.path.join(self.savepath, self.formula_w, stackmode_, cord_)):
                        os.makedirs(os.path.join(self.savepath, self.formula_w, stackmode_, cord_))
                    elif not self.overwrite:
                        mode_count = 2
                        while os.path.exists(os.path.join(self.savepath, self.formula_w, stackmode_add_count, cord_)):
                            stackmode_add_count = f"{stackmode_}_{mode_count}"
                            mode_count += 1
                        os.makedirs(os.path.join(self.savepath, self.formula_w, stackmode_add_count, cord_))
                    savefile = os.path.join(self.savepath, self.formula_w, stackmode_add_count, cord_, f"POSCAR")

                else:
                    if not os.path.exists(self.savepath):
                        os.makedirs(self.savepath)
                    savefile = os.path.join(self.savepath,
                                            f"{self.formula_w}_{stackmode_}_{cord_}-{self.mismatch:.2f}%-{self.natom1}-POSCAR")
                    if not self.overwrite:
                        mode_count = 2
                        while os.path.exists(savefile):
                            savefile = os.path.join(self.savepath,
                                                    f"{self.formula_w}_{stackmode_}{mode_count}_{cord_}-{self.mismatch:.2f}%-{self.natom1}-POSCAR")
                            mode_count += 1

                with open(savefile, 'w') as f:
                    f.write(f"{self.formula_w}-{cord_}\n")
                    f.write("1.0\n")
                    f.write(f"{self.lattice_w[0][0]:.6f} {self.lattice_w[0][1]:.6f} {self.lattice_w[0][2]:.6f}\n")
                    f.write(f"{self.lattice_w[1][0]:.6f} {self.lattice_w[1][1]:.6f} {self.lattice_w[1][2]:.6f}\n")
                    f.write(f"{self.lattice_w[2][0]:.6f} {self.lattice_w[2][1]:.6f} {self.lattice_w[2][2]:.6f}\n")
                    for elem_ in element_pos_dict.keys():
                        f.write(f"{elem_} ")
                    f.write("\n")
                    for elem_ in element_pos_dict.keys():
                        f.write(f"{element_pos_dict[elem_]['count']} ")
                    f.write("\n")
                    f.write("Direct\n")
                    for elem_ in element_pos_dict.keys():
                        for pos_ in element_pos_dict[elem_]['pos']:
                            f.write(f"{pos_[0]:.6f} {pos_[1]:.6f} {pos_[2]:.6f}\n")


class TMDHBilayer:
    def __init__(self, st1, st2=None, a1=None, a2=None, d_intra1=None, d_intra2=None, d_inter=3.3, lv=30.0,
                 tmdtype="2H", savepath=None, savenamemode=1, mismatch_threshold=5.0):
        """
        Constrcting bilayer transition metal dichalcogendies/halides/oxides (TMDs/MHs/TMOs) structures.
        Only hexagonal lattice is supported. The layer group is p3m1, p-3m1 or p-6m2.

        st1: list or str
            The atoms in the first layer or the POSCAR file of the first layer.
            If st1 is a list, for example ['Se', 'Mo', 'S'] for MoSSe.
            Note that the sequence of atoms in the list is ordered by height in Z direction.
            If st1 is a str, it is the path to the POSCAR file of the first layer, and a1, d_intra1 will be ignored.
        st2: list or str
            The atoms in the second layer or the POSCAR file of the second layer.
            If not provided, homogeneous bilayer will be constructed.
        a1: float
            The in-plane lattice constant of the TMDs/MHs/TMOs (Hexagonal lattice) for the first layer.
        a2: float
            The in-plane lattice constant of the TMDs/MHs/TMOs (Hexagonal lattice) for the second layer.
            If not provided, it will be set to a1.
        d_intra1: list
            The distance between two atoms in the same layer for the first layer,
            for example, [1.528, 1.704].
        d_intra2: list
            The distance between two atoms in the same layer for the second layer.
            If not provided, it will be set to d_intra1.
        d_inter: float
            The distance between two atoms in different layers.
        lv: float
            The length of the vacuum layer.
        tmdtype: str
            The type of the TMDs/MHs/TMOs, 1T, 2H, TH (mixed 1T and 2H) or HT (mixed 2H and 1T).
            If st1 is a str, this will be determined by the POSCAR file, and ignore this input and the type of the second layer.
        savepath: str
            The path to save the POSCAR files.
            If not provided, will be saved in the BiPOSCAR_dir directory.
        savenamemode: int
            The mode of the POSCAR file name.
            1: save the POSCAR files in the savepath/formula/genmode/cord*/POSCAR directory.
            2: save the POSCAR files in the savepath directory.
        mismatch_threshold: float
            The maximum tolerance of mismatch_threshold lattice constant between two layers.
        """

        assert tmdtype.upper() in ["1T", "2H", "TH",
                                   "HT"], "Please provide the correct TMD/MH/TMO type: 1T, 2H, TH or HT."

        if savepath is None:
            self.savepath = 'BiPOSCAR_dir'
        else:
            self.savepath = savepath
        if not os.path.exists(self.savepath):
            os.makedirs(self.savepath)
        if savenamemode not in [1, 2]:
            warnings.warn("Please provide the valid savenamemode: 1 or 2. The default is 1.")
            savenamemode = 1
        self.savenamemode = savenamemode
        self.mismatch_threshold = mismatch_threshold
        self.mismatch = 0.0

        self.natom1 = self.ini_atoms()
        lattice_constant = []
        if isinstance(st1, str) or isinstance(st1, Monolayer):
            if isinstance(st1, str):
                mono1 = Monolayer(st1)
            else:
                mono1 = st1
            lattice_constant.append(mono1.checkTMDH(natoms=self.natom1))
            elem1 = mono1.elements
            d_a1 = mono1.get_intra_distances()
            tmd_type1 = mono1.tmdtype
        elif isinstance(st1, list) and all(isinstance(i, str) for i in st1) and len(st1) == self.natom1:
            assert a1 is not None, "Please provide the lattice constant of the TMDs/MHs/TMOs."
            assert d_intra1 is not None and len(
                d_intra1) == self.natom1 - 1, "Please provide the intra distances of the first layer."
            elem1 = st1
            lattice_constant.append(a1)
            d_a1 = d_intra1
            tmd_type1 = tmdtype.upper()
        else:
            raise ValueError("layer1 should be a list of 3 elements or the POSCAR file of the first layer.")

        if st2 is not None:
            if isinstance(st2, str) or isinstance(st2, Monolayer):
                if isinstance(st2, str):
                    mono2 = Monolayer(st2)
                else:
                    mono2 = st2
                lattice_constant.append(mono2.checkTMDH(natoms=self.natom1))
                elem2 = mono2.elements
                d_a2 = mono2.get_intra_distances()
                tmd_type2 = mono2.tmdtype
            elif isinstance(st2, list) and all(isinstance(i, str) for i in st2) and len(st2) == self.natom1:
                assert d_intra2 is not None and len(
                    d_intra2) == self.natom1 - 1, "Please provide the intra distances of the second layer."
                elem2 = st2
                if a2 is not None:
                    lattice_constant.append(a2)
                d_a2 = d_intra2
                tmd_type2 = tmdtype.upper()
            else:
                raise ValueError("layer2 should be a list of 3 elements or the POSCAR file of the second layer.")
        else:
            elem2 = elem1
            lattice_constant.append(lattice_constant[0])
            d_a2 = d_a1
            tmd_type2 = tmd_type1

        # self.elements = elem1 + elem2
        self.check_lattice_mismatch(lattice_constant)
        self.la = np.mean(lattice_constant)
        self.lv = lv
        if tmd_type1 == tmd_type2:
            self.tmdtype = tmd_type1
        else:
            if tmd_type1 == "1T":
                self.tmdtype = "TH"
            elif tmd_type1 == "2H":
                self.tmdtype = "HT"
            else:
                self.tmdtype = tmdtype.upper()

        self.d_intra1 = d_a1
        self.d_intra2 = d_a2
        self.d_inter = d_inter

        e1 = np.array([1, 0])
        e2 = np.array([np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3)])
        self.a_v = e1 * self.la
        self.b_v = e2 * self.la
        self.c_v = (e1 - e2) * self.la / 3
        self.lattice_w = np.array([[self.a_v[0], self.a_v[1], 0], [self.b_v[0], self.b_v[1], 0], [0, 0, self.lv]])

        if elem1 == elem2 or elem1 == elem2[::-1]:
            self.hetero = False
            self.elements = elem1 + elem1
        else:
            self.hetero = True
            self.elements = elem1 + elem2

        if self.tmdtype in ["HT", "TH"]:
            self.hetero = True

        self.set_formula_w()
        if self.savenamemode == 1:
            self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + str(self.natom1)
        self.ini_positions()

    def ini_atoms(self):
        return 3

    def set_formula_w(self):
        if self.elements[0] != self.elements[2] and self.elements[3] != self.elements[-1]:
            self.seqtype = 4
            self.formula_w = self.elements[1] + self.elements[0] + self.elements[2] + self.elements[4] + self.elements[
                3] + self.elements[-1]
        elif self.elements[0] == self.elements[2] and self.elements[3] == self.elements[-1]:
            self.seqtype = 1
            self.formula_w = self.elements[1] + self.elements[0] + str(2) + self.elements[4] + self.elements[3] + str(2)
        elif self.elements[0] != self.elements[2] and self.elements[3] == self.elements[-1]:
            self.seqtype = 2
            self.formula_w = self.elements[1] + self.elements[0] + self.elements[2] + self.elements[4] + self.elements[
                3] + str(2)
        else:
            self.seqtype = 3
            self.formula_w = self.elements[1] + self.elements[0] + str(2) + self.elements[4] + self.elements[3] + \
                             self.elements[-1]

    def ini_positions(self):
        if self.seqtype == 1:
            self.seq_cord = ["cord1"]
        elif self.seqtype == 2:
            self.seq_cord = ["cord1", "cord2"]
        elif self.seqtype == 3:
            self.seq_cord = ["cord1", "cord3"]
        elif self.hetero is False and self.seqtype == 4:
            self.seq_cord = ["cord1", "cord2", "cord3"]
        else:
            self.seq_cord = ["cord1", "cord2", "cord3", "cord4"]

        self.pos_z = {
            "cord1": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[0],
                      self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2],
            "cord2": [self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[1],
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2],
            "cord3": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[0],
                      self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 - self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-2],
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2)],
            "cord4": [self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[1],
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-2],
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2)],
        }

        self.inplane_pos = {"2H": {
            "AA": [2 * self.c_v, self.c_v, 2 * self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v],
            "AA2": [2 * self.c_v, self.c_v, 2 * self.c_v, self.c_v, [0, 0], self.c_v],
            "AA2-2": [self.c_v, [0, 0], self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v],
            "AB": [2 * self.c_v, self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v, self.c_v],
            "AB2": [2 * self.c_v, self.c_v, 2 * self.c_v, 2 * self.c_v, self.a_v, 2 * self.c_v],
            "AB3": [2 * self.c_v, self.c_v, 2 * self.c_v, self.a_v, self.c_v, self.a_v]
        },
            "1T": {
                "AA": [self.c_v, [0, 0], 2 * self.c_v, 2 * self.c_v, [0, 0], self.c_v],
                "AA2": [self.c_v, [0, 0], 2 * self.c_v, self.c_v, 2 * self.c_v, [0, 0]],
                "AA2-2": [2 * self.c_v, [0, 0], self.c_v, [0, 0], 2 * self.c_v, self.c_v],
                "AB": [self.c_v, [0, 0], 2 * self.c_v, self.c_v, [0, 0], 2 * self.c_v],
                "AB2": [2 * self.c_v, self.c_v, self.a_v + self.b_v, self.c_v, [0, 0], 2 * self.c_v],
                "AB3": [self.c_v, [0, 0], 2 * self.c_v, 2 * self.c_v, self.c_v, self.a_v + self.b_v],
            },
            "TH": {
                "AA": [self.c_v, [0, 0], 2 * self.c_v, self.c_v, [0, 0], self.c_v],
                "AA2": [self.c_v, [0, 0], 2 * self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v],
                "AA2-2": [self.c_v, [0, 0], 2 * self.c_v, self.a_v + self.b_v, 2 * self.c_v, self.a_v + self.b_v],
                "AB": [self.c_v, [0, 0], 2 * self.c_v, 2 * self.c_v, [0, 0], 2 * self.c_v],
                "AB2": [self.c_v, [0, 0], 2 * self.c_v, [0, 0], self.c_v, [0, 0]],
                "AB3": [self.c_v, [0, 0], 2 * self.c_v, self.c_v, 2 * self.c_v, self.c_v]
            },
            "HT": {
                "AA": [self.c_v, [0, 0], self.c_v, self.c_v, [0, 0], 2 * self.c_v],
                "AA2": [2 * self.c_v, self.c_v, 2 * self.c_v, self.c_v, [0, 0], 2 * self.c_v],
                "AA2-2": [self.a_v + self.b_v, 2 * self.c_v, self.a_v + self.b_v, self.c_v, [0, 0], 2 * self.c_v],
                "AB": [2 * self.c_v, [0, 0], 2 * self.c_v, self.c_v, [0, 0], 2 * self.c_v],
                "AB2": [[0, 0], self.c_v, [0, 0], self.c_v, [0, 0], 2 * self.c_v],
                "AB3": [self.c_v, 2 * self.c_v, self.c_v, self.c_v, [0, 0], 2 * self.c_v],
            }
        }

    def check_lattice_mismatch(self, la):
        if len(set(la)) > 1:
            self.mismatch = 2 * abs(la[0] - la[1]) / (la[0] + la[1]) * 100
            if self.mismatch > self.mismatch_threshold:
                warnings.warn(f"The lattice mismatch_threshold between the two layers is {self.mismatch:.2f}%, "
                              f"which is larger than the threshold {self.mismatch_threshold}%. \n"
                              f"Please provide a more suitable lattice constant for the two layers.")

    def WritePOSCARsingle(self, stackmode):
        """
        Stacking the bilayer in different stacking modes for 2H or 1T TMDs/MHs/TMOs, and write the POSCAR files.

        genmode: str
            The stacking mode of the bilayer, AA, AA2, AB, AB2, AB3, which correspond to AA, AA', AB, A'B, AB' stacking.
        """
        stackmode = stackmode.upper()
        assert stackmode in ['AA', 'AA2', 'AB', 'AB2', 'AB3', "AA2-2"], "Please provide a valid stacking pattern.\n" \
                                                                        "Valid options are 'AA', 'AA2', 'AB', 'AB2', 'AB3', 'AA2-2'."

        for cord_ in self.seq_cord:
            self.pos = []
            for i, ab in enumerate(self.inplane_pos[self.tmdtype][stackmode]):
                self.pos.append(np.append(ab, self.pos_z[cord_][i]))

            element_pos_dict = {}
            for i, e in enumerate(self.elements):
                if e not in element_pos_dict:
                    element_pos_dict[e] = {}
                    element_pos_dict[e]['count'] = 0
                    element_pos_dict[e]['pos'] = []
                element_pos_dict[e]['count'] += 1
                element_pos_dict[e]['pos'].append(self.pos[i])

            if self.savenamemode == 1:
                if not os.path.exists(os.path.join(self.savepath, self.formula_w, stackmode, cord_)):
                    os.makedirs(os.path.join(self.savepath, self.formula_w, stackmode, cord_))
                savefile = os.path.join(self.savepath, self.formula_w, stackmode, cord_, f"POSCAR")
            else:
                if not os.path.exists(self.savepath):
                    os.makedirs(self.savepath)
                savefile = os.path.join(self.savepath,
                                        f"{self.formula_w}_{stackmode}_{cord_}-{self.mismatch:.2f}%-{str(self.natom1)}-POSCAR")

            with open(savefile, 'w') as f:
                f.write(f"{self.formula_w}-{stackmode}-{cord_}\n")
                f.write("1.0\n")
                f.write(f"{self.lattice_w[0][0]:.6f} {self.lattice_w[0][1]:.6f} {self.lattice_w[0][2]:.6f}\n")
                f.write(f"{self.lattice_w[1][0]:.6f} {self.lattice_w[1][1]:.6f} {self.lattice_w[1][2]:.6f}\n")
                f.write(f"{self.lattice_w[2][0]:.6f} {self.lattice_w[2][1]:.6f} {self.lattice_w[2][2]:.6f}\n")
                for elem_ in element_pos_dict.keys():
                    f.write(f"{elem_} ")
                f.write("\n")
                for elem_ in element_pos_dict.keys():
                    f.write(f"{element_pos_dict[elem_]['count']} ")
                f.write("\n")
                f.write("Cartesian\n")
                for elem_ in element_pos_dict.keys():
                    for pos_ in element_pos_dict[elem_]['pos']:
                        f.write(f"{pos_[0]:.6f} {pos_[1]:.6f} {pos_[2]:.6f}\n")

    def WritePOSCAR(self):
        if self.hetero:
            for stackmode in ['AA', 'AA2', 'AA2-2', 'AB', 'AB2', 'AB3']:
                self.WritePOSCARsingle(stackmode)
        else:
            for stackmode in ['AA', 'AA2', 'AB', 'AB2', 'AB3']:
                self.WritePOSCARsingle(stackmode)


class MNXYBilayer(TMDHBilayer):
    """
    Construction of bilayers for the structure of Al2S2 whose layer group is P-3m1(1T) or P-6m2 (2H).
    Only hexagonal lattices are supported.
    """

    def ini_atoms(self):
        return 4

    def set_formula_w(self):
        if (self.elements[0] != self.elements[3] or self.elements[1] != self.elements[2]) and \
                (self.elements[4] != self.elements[-1] or self.elements[5] != self.elements[-2]):
            self.seqtype = 4
            self.formula_w = self.elements[1] + self.elements[2] + self.elements[0] + self.elements[3] + \
                             self.elements[5] + self.elements[6] + self.elements[4] + self.elements[-1]
        elif self.elements[0] == self.elements[3] and self.elements[1] == self.elements[2] and self.elements[4] == \
                self.elements[-1] and self.elements[5] == self.elements[-2]:
            self.seqtype = 1
            self.formula_w = self.elements[1] + str(2) + self.elements[0] + str(2) + self.elements[5] + str(2) \
                             + self.elements[4] + str(2)
        elif (self.elements[0] != self.elements[3] or self.elements[1] != self.elements[2]) and self.elements[4] == \
                self.elements[-1] and self.elements[5] == self.elements[-2]:
            self.seqtype = 2
            self.formula_w = self.elements[1] + self.elements[2] + self.elements[0] + self.elements[3] \
                             + self.elements[5] + str(2) + self.elements[4] + str(2)
        else:
            self.seqtype = 3
            self.formula_w = self.elements[1] + str(2) + self.elements[0] + str(2) + self.elements[5] + \
                             self.elements[6] + self.elements[4] + self.elements[-1]

    def ini_positions(self):
        if self.seqtype == 1:
            self.seq_cord = ["cord1"]
        elif self.seqtype == 2:
            self.seq_cord = ["cord1", "cord2"]
        elif self.seqtype == 3:
            self.seq_cord = ["cord1", "cord3"]
        elif self.hetero is False and self.seqtype == 4:
            self.seq_cord = ["cord1", "cord2", "cord3"]
        else:
            self.seq_cord = ["cord1", "cord2", "cord3", "cord4"]

        self.pos_z = {
            "cord1": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[0],
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[0] + self.d_intra1[1],
                      self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1] - self.d_intra2[-2],
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2],
            "cord2": [self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[1] + self.d_intra1[2],
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[2],
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1] - self.d_intra2[-2],
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2],
            "cord3": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[0],
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[0] + self.d_intra1[1],
                      self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 - self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[0],
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[0] - self.d_intra2[1],
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2)],
            "cord4": [self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[1] + self.d_intra1[2],
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1[2],
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[0],
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[0] - self.d_intra2[1],
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2)],
        }

        self.inplane_pos = {"2H": {
            "AA": [2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v, 2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v],
            "AA2": [2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v, self.c_v, [0, 0], [0, 0], self.c_v],
            "AA2-2": [self.c_v, [0, 0], [0, 0], self.c_v, 2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v],
            "AB": [2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v, 2 * self.c_v, self.c_v],
            "AB2": [2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v, 2 * self.c_v, self.a_v, self.a_v, 2 * self.c_v],
            "AB3": [2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v, self.a_v, self.c_v, self.c_v, self.a_v]
        },
            "1T": {
                "AA": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, 2 * self.c_v, [0, 0], [0, 0], self.c_v],
                "AA2": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, self.c_v, 2 * self.c_v, 2 * self.c_v, [0, 0]],
                "AA2-2": [2 * self.c_v, [0, 0], [0, 0], self.c_v, [0, 0], 2 * self.c_v, 2 * self.c_v, self.c_v],
                "AB": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, self.c_v, [0, 0], [0, 0], 2 * self.c_v],
                "AB2": [2 * self.c_v, self.c_v, self.c_v, self.a_v + self.b_v, self.c_v, [0, 0], [0, 0], 2 * self.c_v],
                "AB3": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, 2 * self.c_v, self.c_v, self.c_v, self.a_v + self.b_v],
            },
            "TH": {
                "AA": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, self.c_v, [0, 0], [0, 0], self.c_v],
                "AA2": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, 2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v],
                "AA2-2": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, self.a_v + self.b_v, 2 * self.c_v, 2 * self.c_v,
                          self.a_v + self.b_v],
                "AB": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, 2 * self.c_v, [0, 0], [0, 0], 2 * self.c_v],
                "AB2": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, [0, 0], self.c_v, self.c_v, [0, 0]],
                "AB3": [self.c_v, [0, 0], [0, 0], 2 * self.c_v, self.c_v, 2 * self.c_v, 2 * self.c_v, self.c_v]
            },
            "HT": {
                "AA": [self.c_v, [0, 0], [0, 0], self.c_v, self.c_v, [0, 0], [0, 0], 2 * self.c_v],
                "AA2": [2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v, self.c_v, [0, 0], [0, 0], 2 * self.c_v],
                "AA2-2": [self.a_v + self.b_v, 2 * self.c_v, 2 * self.c_v, self.a_v + self.b_v, self.c_v, [0, 0],
                          [0, 0], 2 * self.c_v],
                "AB": [2 * self.c_v, [0, 0], [0, 0], 2 * self.c_v, self.c_v, [0, 0], [0, 0], 2 * self.c_v],
                "AB2": [[0, 0], self.c_v, self.c_v, [0, 0], self.c_v, [0, 0], [0, 0], 2 * self.c_v],
                "AB3": [self.c_v, 2 * self.c_v, 2 * self.c_v, self.c_v, self.c_v, [0, 0], [0, 0], 2 * self.c_v],
            }
        }


class N2Bilayer:
    def __init__(self, st1, st2=None, la1=None, la2=None, d_intra1=None, d_intra2=None, d_inter=3.3, lv=25.0,
                 savepath=None, savenamemode=1, mismatch_threshold=5.0):
        """
        Construct the N2 bilayer with the given layers and positions, where the N2 means the layer is composed of two atoms.
        Currently, only the hexagonal lattice is supported.

        st1: list or str
            The atoms in the first layer or the POSCAR file of the first layer.
            If st1 is a list, note that the sequence of atoms in the list is ordered by height in Z direction.
            If st1 is a str, it is the path to the POSCAR file of the first layer, and a1, d_intra1 will be ignored.
        st2: list or str
            The atoms in the second layer or the POSCAR file of the second layer.
            If not provided, homogeneous bilayer will be constructed.
        la1: float
            The in-plane lattice constant for the first layer.
        la2: float
            The in-plane lattice constant for the second layer.
            If not provided, it will be set to a1.
        d_intra1: list
            The distance between two atoms in the first layer.
        d_intra2: list
            The distance between two atoms in the second layer.
            If not provided, it will be set to d_intra1.
        d_inter: float
            The distance between two atoms in different layers.
        lv: float
            The length of the vacuum layer.
        savepath: str
            The path to save the POSCAR files.
            If not provided, will be saved in the BiPOSCAR_dir directory.
        savenamemode: int
            The mode to name the saved POSCAR files.
            1: The POSCAR files are saved in the folder f"{savepath}/{formula}_{mismatch}/{genmode}/{cord*}/POSCAR".
            2: All the POSCAR files are saved in the same folder f"{savepath}/{formula}_{genmode}_{cord*}-{mismatch}-POSCAR".
        mismatch_threshold: float
            Maximum tolerance of mismatch_threshold lattice constant between two layers.
        """
        if savepath is None:
            self.savepath = 'BiPOSCAR_dir'
        else:
            self.savepath = savepath
        if not os.path.exists(self.savepath):
            os.makedirs(self.savepath)
        self.savenamemode = savenamemode

        self.mismatch_threshold = mismatch_threshold
        self.mismatch = 0.0

        self.tmdtype = "2H"

        self.natom1 = 2
        elem1, Latype1, a1, d_a1, elem2, Latype2, a2, d_a2 = self.ini_st(st1, st2, la1, la2, d_intra1, d_intra2)

        assert Latype1 == Latype2 and Latype1 == "Hexagonal", "Currently, only the hexagonal lattice is supported."
        if self.hetero:
            self.mismatch = np.abs(a1 - a2) / a1 * 100
            if self.mismatch > self.mismatch_threshold:
                warnings.warn(
                    f"The mismatch between two layers is {self.mismatch:.2f}%, which is larger than the threshold "
                    f"{self.mismatch_threshold}%. Please check the input parameters.")
            self.la = (a1 + a2) / 2
            self.elements = elem1 + elem2
        else:
            self.la = a1
            self.elements = elem1 + elem1
        self.lv = lv
        self.d_intra1 = d_a1
        self.d_intra2 = d_a2
        self.d_inter = d_inter

        e1 = np.array([1, 0])
        e2 = np.array([np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3)])
        self.a_v = e1 * self.la
        self.b_v = e2 * self.la
        self.c_v = (e1 - e2) * self.la / 3
        self.lattice_w = np.array([[self.a_v[0], self.a_v[1], 0], [self.b_v[0], self.b_v[1], 0], [0, 0, self.lv]])

        self.set_formula_w()
        self.ini_positions()

    def ini_st(self, st1, st2, la1, la2, d_intra1, d_intra2):
        if isinstance(st1, str) or isinstance(st1, Monolayer):
            if isinstance(st1, str):
                mono1 = Monolayer(st1)
            else:
                mono1 = st1
            Latype1 = mono1.checkLatticeType()
            a1 = mono1.a
            elem1 = mono1.elements
            d_a1 = np.diff(mono1.positions[:, 2])
        elif isinstance(st1, list):
            assert len(st1) == self.natom1, "The first layer should be composed of two atoms."
            assert isinstance(la1, float), "Please provide the in-plane lattice constant for the first layer."
            assert isinstance(d_intra1, list) and len(
                d_intra1) == self.natom1 - 1, "Please provide the distance between two atoms in the first layer."
            elem1 = st1
            Latype1, a1 = "Hexagonal", la1
            d_a1 = d_intra1
        else:
            raise ValueError("Please provide the atoms in the first layer or the POSCAR file of the first layer.")

        self.hetero = False
        if (isinstance(st2, str) or isinstance(st2, Monolayer)) and st1 != st2:
            if isinstance(st2, str):
                mono2 = Monolayer(st2)
            else:
                mono2 = st2
            elem2 = mono2.elements
            if elem1 == elem2 or elem1 == elem2[::-1]:
                elem2, Latype2, a2, d_a2 = elem1, Latype1, a1, d_a1
            else:
                self.hetero = True
                Latype2 = mono2.checkLatticeType()
                a2 = mono2.a
                d_a2 = np.diff(mono2.positions[:, 2])
        elif isinstance(st2, list):
            if st2 == elem1 or st2 == elem1[::-1]:
                elem2, Latype2, a2, d_a2 = elem1, Latype1, a1, d_a1
            else:
                self.hetero = True
                elem2 = st2
                if isinstance(la2, float):
                    Latype2, a2 = "Hexagonal", la2
                else:
                    Latype2, a2 = Latype1, a1

                if isinstance(d_intra2, list) and len(d_intra2) == self.natom1 - 1:
                    d_a2 = d_intra2
                else:
                    d_a2 = d_a1
        else:
            elem2, Latype2, a2, d_a2 = elem1, Latype1, a1, d_a1
        return elem1, Latype1, a1, d_a1, elem2, Latype2, a2, d_a2

    def set_formula_w(self):
        if self.elements[0] != self.elements[1] and self.elements[2] != self.elements[3]:
            self.seqtype = 4
            self.formula_w = "".join(self.elements)
        elif self.elements[0] == self.elements[1] and self.elements[2] == self.elements[3]:
            self.seqtype = 1
            self.formula_w = self.elements[0] + str(2) + self.elements[2] + str(2)
        elif self.elements[0] != self.elements[1] and self.elements[2] == self.elements[3]:
            self.seqtype = 2
            self.formula_w = self.elements[0] + self.elements[1] + self.elements[2] + str(2)
        else:
            self.seqtype = 3
            self.formula_w = self.elements[0] + str(2) + self.elements[2] + self.elements[3]
        if self.savenamemode == 1:
            self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + str(self.natom1)

    def ini_positions(self):
        self.inplane_pos = {
            "AA": [self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v],
            "AB": [[0, 0], self.c_v, self.c_v, 2 * self.c_v]
        }

        self.outplane_pos = {
            "cord1": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2,
                      self.lv / 2 - self.d_inter / 2],
            "cord2": [self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2,
                      self.lv / 2 - self.d_inter / 2],
            "cord3": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 - self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2],
            "cord4": [self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2]
        }

        def pos_gen(stacktype, cordtype):
            pos_ = []
            for i in range(4):
                pos_.append(np.append(self.inplane_pos[stacktype][i], self.outplane_pos[cordtype][i]))
            return pos_

        self.stackpos = {"AA": pos_gen("AA", "cord1"),
                         "AA2-1": pos_gen("AA", "cord2"),
                         "AA2-2": pos_gen("AA", "cord3"),
                         "AB": pos_gen("AB", "cord1"),
                         "AB2-1": pos_gen("AB", "cord2"),
                         "AB2-2": pos_gen("AB", "cord3"),
                         "AB3": pos_gen("AB", "cord4")}

        homo_delta_h = True
        if abs(self.d_intra1 - 0.0) < 1e-3 and abs(self.d_intra2 - 0.0) < 1e-3:
            self.stackmodes = ["AA", "AB"]
            homo_delta_h = False
        elif abs(self.d_intra1 - 0.0) >= 1e-3 and abs(self.d_intra2 - 0.0) < 1e-3:
            self.stackmodes = ["AA", "AA2-1", "AB", "AB2-1"]
        elif abs(self.d_intra1 - 0.0) < 1e-3 and abs(self.d_intra2 - 0.0) >= 1e-3:
            self.stackmodes = ["AA", "AA2-2", "AB", "AB2-2"]
        else:
            self.stackmodes = ["AA", "AA2-1", "AB", "AB2-1", "AB3"]

        self.stackelements = {}
        if self.seqtype == 1:
            elem0 = self.elements[:]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0}
        elif self.seqtype == 2:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[1] = elem1[1], elem1[0]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1}
            self.stackelements["AA"] = {"cord1": elem0}
        elif self.seqtype == 3:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[2], elem1[3] = elem1[3], elem1[2]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1}
            self.stackelements["AA"] = {"cord1": elem0}
        else:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[1] = elem1[1], elem1[0]
            elem2 = elem0[:]
            elem2[2], elem2[3] = elem2[3], elem2[2]
            elem3 = elem0[:]
            elem3[0], elem3[1], elem3[2], elem3[3] = elem3[1], elem3[0], elem3[3], elem3[2]
            if self.hetero:
                for stack_mode in self.stackmodes:
                    self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1, "cord3": elem2, "cord4": elem3}
            else:
                if homo_delta_h:
                    self.stackelements = {"AA": {"cord1": elem0, "cord2": elem1, "cord3": elem2},
                                          "AA2-1": {"cord1": elem0, "cord2": elem1, "cord3": elem3},
                                          "AB": {"cord1": elem0, "cord2": elem1, "cord3": elem2},
                                          "AB2-1": {"cord1": elem0, "cord2": elem1, "cord3": elem2},
                                          "AB3": {"cord1": elem0, "cord2": elem1, "cord3": elem2}}
                else:
                    self.stackelements = {"AA": {"cord1": elem0, "cord2": elem1},
                                          "AA2-1": {"cord1": elem0, "cord2": elem1, "cord3": elem3},
                                          "AB": {"cord1": elem0, "cord2": elem1, "cord3": elem2},
                                          "AB2-1": {"cord1": elem0, "cord2": elem1, "cord3": elem2},
                                          "AB3": {"cord1": elem0, "cord2": elem1, "cord3": elem2}}

    def poszc(self, cord_, pos_):
        return pos_

    def WritePOSCAR(self):
        """
        Write the POSCAR file for the given stacking pattern.
        """
        for stackmode in self.stackmodes:
            for cord_ in self.stackelements[stackmode]:
                element_pos_dict = {}
                pos_ = self.stackpos[stackmode]
                pos_ = self.poszc(cord_, pos_)
                for i, e in enumerate(self.stackelements[stackmode][cord_]):
                    if e not in element_pos_dict:
                        element_pos_dict[e] = {}
                        element_pos_dict[e]['count'] = 0
                        element_pos_dict[e]['pos'] = []
                    element_pos_dict[e]['count'] += 1
                    element_pos_dict[e]['pos'].append(pos_[i])

                element_pos_dict = dict(sorted(element_pos_dict.items(), key=lambda x: self.elements.index(x[0])))

                if self.savenamemode == 1:
                    if not os.path.exists(os.path.join(self.savepath, self.formula_w, stackmode, cord_)):
                        os.makedirs(os.path.join(self.savepath, self.formula_w, stackmode, cord_))
                    savefile = os.path.join(self.savepath, self.formula_w, stackmode, cord_, f"POSCAR")
                else:
                    if not os.path.exists(self.savepath):
                        os.makedirs(self.savepath)
                    savefile = os.path.join(self.savepath,
                                            f"{self.formula_w}_{stackmode}_{cord_}-{self.mismatch:.2f}%-{str(self.natom1)}-POSCAR")

                with open(savefile, 'w') as f:
                    f.write(f"{self.formula_w}-{stackmode}-{cord_}\n")
                    f.write("1.0\n")
                    f.write(f"{self.lattice_w[0][0]:.6f} {self.lattice_w[0][1]:.6f} {self.lattice_w[0][2]:.6f}\n")
                    f.write(f"{self.lattice_w[1][0]:.6f} {self.lattice_w[1][1]:.6f} {self.lattice_w[1][2]:.6f}\n")
                    f.write(f"{self.lattice_w[2][0]:.6f} {self.lattice_w[2][1]:.6f} {self.lattice_w[2][2]:.6f}\n")
                    for elem_ in element_pos_dict.keys():
                        f.write(f"{elem_} ")
                    f.write("\n")
                    for elem_ in element_pos_dict.keys():
                        f.write(f"{element_pos_dict[elem_]['count']} ")
                    f.write("\n")
                    f.write("Cartesian\n")
                    for elem_ in element_pos_dict.keys():
                        for pos_ in element_pos_dict[elem_]['pos']:
                            f.write(f"{pos_[0]:.6f} {pos_[1]:.6f} {pos_[2]:.6f}\n")


class N2TMDHbilayer(N2Bilayer):
    """
    Construct a bilayer system with the one layer being a N2 monolaer and the other layer being a TMD/MH/TMO monolayer.
    For example, Graphene-MoS2, etc.
    Only support the 2H-TMDs/MHs/TMOs currently.
    """

    def ini_st(self, st1, st2, la1, la2, d_intra1, d_intra2):
        layer_2atom = 1

        if isinstance(st1, str) or isinstance(st1, Monolayer):
            if isinstance(st1, str):
                mono1 = Monolayer(st1)
            else:
                mono1 = st1
            Latype1 = mono1.checkLatticeType()
            a1 = mono1.a
            elem1 = mono1.elements
            d_a1 = np.diff(mono1.positions[:, 2])
            if len(elem1) == 3:
                layer_2atom = 2
                a_ = mono1.checkTMDH()
                self.tmdtype = mono1.tmdtype
        elif isinstance(st1, list):
            assert isinstance(la1, float), "Please provide the in-plane lattice constant for the first layer."
            assert isinstance(d_intra1, list), "Please provide the interlayer distance for the first layer."
            elem1 = st1
            Latype1, a1 = "Hexagonal", la1
            d_a1 = d_intra1
            if len(elem1) == 3:
                layer_2atom = 2
                warnings.warn(
                    "The TMD/MH/TMO type is set to 2H by default. If not, please provide the POSCAR file of the first layer.")
        else:
            raise ValueError("Please provide the atoms in the first layer or the POSCAR file of the first layer.")

        if isinstance(st2, str) or isinstance(st2, Monolayer):
            if isinstance(st2, str):
                mono2 = Monolayer(st2)
            else:
                mono2 = st2
            Latype2 = mono2.checkLatticeType()
            a2 = mono2.a
            elem2 = mono2.elements
            d_a2 = np.diff(mono2.positions[:, 2])
            if len(elem2) == 3:
                a_ = mono2.checkTMDH()
                self.tmdtype = mono2.tmdtype
                if len(elem1) != 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD/MH/TMO monolayer.")
            elif len(elem2) == 2:
                if len(elem1) == 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD/MH/TMO monolayer.")
            else:
                raise ValueError("Please provide the correct atoms in the second layer.")
        elif isinstance(st2, list):
            assert isinstance(la2, float), "Please provide the in-plane lattice constant for the first layer."
            assert isinstance(d_intra2, list), "Please provide the interlayer distance for the first layer."
            elem2 = st2
            Latype2, a2 = "Hexagonal", la2
            d_a2 = d_intra2
            if len(elem2) == 3:
                warnings.warn(
                    "The TMD/MH/TMO type is set to 2H by default. If not, please provide the POSCAR file of the first layer.")
                if len(elem1) != 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD monolayer.")
            elif len(elem2) == 2:
                if len(elem1) == 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD monolayer.")
            else:
                raise ValueError("Please provide the correct atoms in the second layer.")
        else:
            raise ValueError("Please provide the atoms in the second layer or the POSCAR file of the first layer.")
        self.hetero = True

        if self.tmdtype == "1T":
            raise ValueError("This class only supports the 2H-TMDs/MHs/TMOs currently.")

        if layer_2atom == 1:
            return elem1, Latype1, a1, d_a1, elem2, Latype2, a2, d_a2
        else:
            return elem2, Latype2, a2, d_a2, elem1, Latype1, a1, d_a1

    def set_formula_w(self):
        if self.elements[0] != self.elements[1] and self.elements[2] != self.elements[4]:
            self.seqtype = 4
            self.formula_w = "".join(self.elements)
        elif self.elements[0] == self.elements[1] and self.elements[2] == self.elements[4]:
            self.seqtype = 1
            self.formula_w = self.elements[0] + str(2) + self.elements[3] + self.elements[2] + str(2)
        elif self.elements[0] != self.elements[1] and self.elements[2] == self.elements[4]:
            self.seqtype = 2
            self.formula_w = self.elements[0] + self.elements[1] + self.elements[3] + self.elements[2] + str(2)
        else:
            self.seqtype = 3
            self.formula_w = self.elements[0] + str(2) + self.elements[3] + self.elements[2] + self.elements[4]
        if self.savenamemode == 1:
            self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + str(self.natom1)

    def ini_positions(self):
        self.inplane_pos = {
            "AA": [self.c_v, 2 * self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v],
            "AB": [[0, 0], self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v],
            "AC": [[0, 0], 2 * self.c_v, 2 * self.c_v, self.c_v, 2 * self.c_v],
        }

        self.outplane_pos = {
            "cord1": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2],
            "cord2": [self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2]
        }

        def pos_gen(stacktype, cordtype):
            pos_ = []
            for i in range(5):
                pos_.append(np.append(self.inplane_pos[stacktype][i], self.outplane_pos[cordtype][i]))
            return pos_

        self.stackpos = {"AA": pos_gen("AA", "cord1"),
                         "AA2": pos_gen("AA", "cord2"),
                         "AB": pos_gen("AB", "cord1"),
                         "AB2": pos_gen("AB", "cord2"),
                         "AC": pos_gen("AC", "cord1"),
                         "AC2": pos_gen("AC", "cord2")}

        if abs(self.d_intra1 - 0.0) >= 1e-3:
            self.stackmodes = ["AA", "AA2", "AB", "AB2", "AC", "AC2"]
        else:
            self.stackmodes = ["AA", "AB", "AC"]

        self.stackelements = {}
        self.poszc_flag = False
        if self.seqtype == 1:
            elem0 = self.elements[:]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0}
        elif self.seqtype == 2:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[1] = elem1[1], elem1[0]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1}
        elif self.seqtype == 3:
            self.poszc_flag = True
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[2], elem1[4] = elem1[4], elem1[2]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord3": elem1}
        else:
            self.poszc_flag = True
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[1] = elem1[1], elem1[0]
            elem2 = elem0[:]
            elem2[2], elem2[4] = elem2[4], elem2[2]
            elem3 = elem0[:]
            elem3[0], elem3[1], elem3[2], elem3[4] = elem3[1], elem3[0], elem3[4], elem3[2]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1, "cord3": elem2, "cord4": elem3}

    def poszc(self, cord_, pos_):
        if cord_ == "cord3" or cord_ == "cord4":
            if self.poszc_flag:
                pos_[-2][2] = self.lv / 2 - self.d_inter / 2 - self.d_intra2[0]
        return pos_


class N2MNXYBilayer(N2Bilayer):
    """
    Construct a bilayer system with the one layer being a N2 monolaer and the other layer being a MNXY monolayer.
    For example, Graphene-Al2S2, etc.
    Only support the 2H-MNXY currently.
    """

    def ini_st(self, st1, st2, la1, la2, d_intra1, d_intra2):
        layer_2atom = 1

        if isinstance(st1, str) or isinstance(st1, Monolayer):
            if isinstance(st1, str):
                mono1 = Monolayer(st1)
            else:
                mono1 = st1
            Latype1 = mono1.checkLatticeType()
            a1 = mono1.a
            elem1 = mono1.elements
            d_a1 = np.diff(mono1.positions[:, 2])
            if len(elem1) == 4:
                layer_2atom = 2
                a_ = mono1.checkTMDH(natoms=4)
                self.tmdtype = mono1.tmdtype
        elif isinstance(st1, list):
            assert isinstance(la1, float), "Please provide the in-plane lattice constant for the first layer."
            assert isinstance(d_intra1, list), "Please provide the interlayer distance for the first layer."
            elem1 = st1
            Latype1, a1 = "Hexagonal", la1
            d_a1 = d_intra1
            if len(elem1) == 4:
                layer_2atom = 2
                warnings.warn(
                    "The TMD/MH/TMO type is set to 2H by default. If not, please provide the POSCAR file of the first layer.")
        else:
            raise ValueError("Please provide the atoms in the first layer or the POSCAR file of the first layer.")

        if isinstance(st2, str) or isinstance(st2, Monolayer):
            if isinstance(st2, str):
                mono2 = Monolayer(st2)
            else:
                mono2 = st2
            Latype2 = mono2.checkLatticeType()
            a2 = mono2.a
            elem2 = mono2.elements
            d_a2 = np.diff(mono2.positions[:, 2])
            if len(elem2) == 4:
                a_ = mono2.checkTMDH(natoms=4)
                self.tmdtype = mono2.tmdtype
                if len(elem1) != 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD/MH/TMO monolayer.")
            elif len(elem2) == 2:
                if len(elem1) == 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD/MH/TMO monolayer.")
            else:
                raise ValueError("Please provide the correct atoms in the second layer.")
        elif isinstance(st2, list):
            assert isinstance(la2, float), "Please provide the in-plane lattice constant for the first layer."
            assert isinstance(d_intra2, list), "Please provide the interlayer distance for the first layer."
            elem2 = st2
            Latype2, a2 = "Hexagonal", la2
            d_a2 = d_intra2
            if len(elem2) == 4:
                warnings.warn(
                    "The TMD/MH/TMO type is set to 2H by default. If not, please provide the POSCAR file of the first layer.")
                if len(elem1) != 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD/MH/TMO monolayer.")
            elif len(elem2) == 2:
                if len(elem1) == 2:
                    raise ValueError(
                        "This class is used to construct a bilayer system with one layer being a N2 monolayer and the other layer being a TMD/MH/TMO monolayer.")
            else:
                raise ValueError("Please provide the correct atoms in the second layer.")
        else:
            raise ValueError("Please provide the atoms in the second layer or the POSCAR file of the first layer.")
        self.hetero = True

        if self.tmdtype == "1T":
            raise ValueError("This class only supports the 2H-MNXY currently.")

        if layer_2atom == 1:
            return elem1, Latype1, a1, d_a1, elem2, Latype2, a2, d_a2
        else:
            return elem2, Latype2, a2, d_a2, elem1, Latype1, a1, d_a1

    def set_formula_w(self):
        if self.elements[0] != self.elements[1] and (
                self.elements[2] != self.elements[5] or self.elements[3] != self.elements[4]):
            self.seqtype = 4
            self.formula_w = "".join(self.elements)
        elif self.elements[0] == self.elements[1] and self.elements[2] == self.elements[5] and self.elements[3] == \
                self.elements[4]:
            self.seqtype = 1
            self.formula_w = self.elements[0] + str(2) + self.elements[3] + str(2) + self.elements[2] + str(2)
        elif self.elements[0] != self.elements[1] and self.elements[2] == self.elements[5] and self.elements[3] == \
                self.elements[4]:
            self.seqtype = 2
            self.formula_w = self.elements[0] + self.elements[1] + self.elements[3] + str(2) + self.elements[2] + str(2)
        else:
            self.seqtype = 3
            self.formula_w = self.elements[0] + str(2) + self.elements[3] + self.elements[4] + self.elements[2] + \
                             self.elements[5]
        if self.savenamemode == 1:
            self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + str(self.natom1)

    def ini_positions(self):
        self.inplane_pos = {
            "AA": [self.c_v, 2 * self.c_v, 2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v],
            "AB": [[0, 0], self.c_v, 2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v],
            "AC": [[0, 0], 2 * self.c_v, 2 * self.c_v, self.c_v, self.c_v, 2 * self.c_v],
        }

        self.outplane_pos = {
            "cord1": [self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1] - self.d_intra2[-2],
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2],
            "cord2": [self.lv / 2 + self.d_inter / 2 + self.d_intra1,
                      self.lv / 2 + self.d_inter / 2,
                      self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2),
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1] - self.d_intra2[-2],
                      self.lv / 2 - self.d_inter / 2 - self.d_intra2[-1],
                      self.lv / 2 - self.d_inter / 2]
        }

        def pos_gen(stacktype, cordtype):
            pos_ = []
            for i in range(6):
                pos_.append(np.append(self.inplane_pos[stacktype][i], self.outplane_pos[cordtype][i]))
            return pos_

        self.stackpos = {"AA": pos_gen("AA", "cord1"),
                         "AA2": pos_gen("AA", "cord2"),
                         "AB": pos_gen("AB", "cord1"),
                         "AB2": pos_gen("AB", "cord2"),
                         "AC": pos_gen("AC", "cord1"),
                         "AC2": pos_gen("AC", "cord2")}

        if abs(self.d_intra1 - 0.0) >= 1e-3:
            self.stackmodes = ["AA", "AA2", "AB", "AB2", "AC", "AC2"]
        else:
            self.stackmodes = ["AA", "AB", "AC"]

        self.stackelements = {}
        self.poszc_flag = False
        if self.seqtype == 1:
            elem0 = self.elements[:]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0}
        elif self.seqtype == 2:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[1] = elem1[1], elem1[0]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1}
        elif self.seqtype == 3:
            self.poszc_flag = True
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[2], elem1[3], elem1[4], elem1[5] = elem1[5], elem1[4], elem1[3], elem1[2]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord3": elem1}
        else:
            self.poszc_flag = True
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[1] = elem1[1], elem1[0]
            elem2 = elem0[:]
            elem2[2], elem2[3], elem2[4], elem2[5] = elem2[5], elem2[4], elem2[3], elem2[2]
            elem3 = elem0[:]
            elem3[0], elem3[1], elem3[2], elem3[3], elem3[4], elem3[5] = elem3[1], elem3[0], elem3[5], elem3[4], elem3[
                3], elem3[2]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1, "cord3": elem2, "cord4": elem3}

    def poszc(self, cord_, pos_):
        if cord_ == "cord3" or cord_ == "cord4":
            if self.poszc_flag:
                pos_[-3][2] = self.lv / 2 - self.d_inter / 2 - self.d_intra2[0] - self.d_intra2[1]
                pos_[-2][2] = self.lv / 2 - self.d_inter / 2 - self.d_intra2[0]
        return pos_


class TMDHBilayerSquare:
    def __init__(self, st1, st2=None, la1=None, la2=None, d_intra1=None, d_intra2=None, d_inter=3.3, lv=25.0,
                 savepath=None, savenamemode=1, mismatch_threshold=5.0):
        """
        Constrcting bilayer transition metal dichalcogendies/halides/oxides (TMDs/MHs/TMOs) structures, whose layer group is p-4m2.
        Only the square lattice is supported.

        st1: list or str
            The atoms in the first layer or the POSCAR file of the first layer.
            If st1 is a list, note that the sequence of atoms in the list is ordered by height in Z direction.
            If st1 is a str, it is the path to the POSCAR file of the first layer, and a1, d_intra1 will be ignored.
        st2: list or str
            The atoms in the second layer or the POSCAR file of the second layer.
            If not provided, homogeneous bilayer will be constructed.
        la1: float
            The in-plane lattice constant for the first layer.
        la2: float
            The in-plane lattice constant for the second layer.
            If not provided, it will be set to a1.
        d_intra1: list
            The distance between two atoms in the first layer.
        d_intra2: list
            The distance between two atoms in the second layer.
            If not provided, it will be set to d_intra1.
        d_inter: float
            The distance between two atoms in different layers.
        lv: float
            The length of the vacuum layer.
        savepath: str
            The path to save the POSCAR files.
            If not provided, will be saved in the BiPOSCAR_dir directory.
        savenamemode: int
            The mode to name the saved POSCAR files.
            1: The POSCAR files are saved in the folder f"{savepath}/{formula}_{mismatch}/{genmode}/{cord*}/POSCAR".
            2: All the POSCAR files are saved in the same folder f"{savepath}/{formula}_{genmode}_{cord*}-{mismatch}-POSCAR".
        mismatch_threshold: float
            Maximum tolerance of mismatch_threshold lattice constant between two layers.
        """
        if savepath is None:
            self.savepath = 'BiPOSCAR_dir'
        else:
            self.savepath = savepath
        if not os.path.exists(self.savepath):
            os.makedirs(self.savepath)
        self.savenamemode = savenamemode

        self.mismatch_threshold = mismatch_threshold
        self.mismatch = 0.0

        self.ini_atoms()
        elem1, Latype1, a1, d_a1, elem2, Latype2, a2, d_a2 = self.ini_st(st1, st2, la1, la2, d_intra1, d_intra2)

        assert Latype1 == Latype2 and Latype1 == "Square", "The lattice type of the two layers should be the same and be Square."
        if self.hetero:
            self.mismatch = np.abs(a1 - a2) / a1 * 100
            if self.mismatch > self.mismatch_threshold:
                warnings.warn(
                    f"The mismatch between two layers is {self.mismatch:.2f}%, which is larger than the threshold "
                    f"{self.mismatch_threshold}%. Please check the input parameters.")
            self.la = (a1 + a2) / 2
            self.elements = elem1 + elem2
        else:
            self.la = a1
            self.elements = elem1 + elem1
        self.lv = lv
        self.d_intra1 = d_a1
        self.d_intra2 = d_a2
        self.d_inter = d_inter

        self.lattice_w = np.array([[self.la, 0, 0], [0, self.la, 0], [0, 0, self.lv]])

        self.set_formula_w()
        self.ini_positions()

    def ini_atoms(self):
        self.natom = 3

    def ini_st(self, st1, st2, la1, la2, d_intra1, d_intra2):
        if isinstance(st1, str) or isinstance(st1, Monolayer):
            if isinstance(st1, str):
                mono1 = Monolayer(st1)
            else:
                mono1 = st1
            Latype1 = mono1.checkLatticeType()
            a1 = mono1.a
            elem1 = mono1.elements
            d_a1 = np.diff(mono1.positions[:, 2])
        elif isinstance(st1, list):
            assert len(st1) == self.natom, "The first layer should be composed of two atoms."
            assert isinstance(la1, float), "Please provide the in-plane lattice constant for the first layer."
            assert isinstance(d_intra1, list) and len(
                d_intra1) == self.natom - 1, "Please provide the distance between two atoms in the first layer."
            elem1 = st1
            Latype1, a1 = "Square", la1
            d_a1 = d_intra1
        else:
            raise ValueError("Please provide the atoms in the first layer or the POSCAR file of the first layer.")

        self.hetero = False
        if (isinstance(st2, str) or isinstance(st2, Monolayer)) and st1 != st2:
            if isinstance(st2, str):
                mono2 = Monolayer(st2)
            else:
                mono2 = st2
            elem2 = mono2.elements
            if elem1 == elem2 or elem1 == elem2[::-1]:
                elem2, Latype2, a2, d_a2 = elem1, Latype1, a1, d_a1
            else:
                self.hetero = True
                Latype2 = mono2.checkLatticeType()
                a2 = mono2.a
                d_a2 = np.diff(mono2.positions[:, 2])
        elif isinstance(st2, list):
            if st2 == elem1 or st2 == elem1[::-1]:
                elem2, Latype2, a2, d_a2 = elem1, Latype1, a1, d_a1
            else:
                self.hetero = True
                elem2 = st2
                if isinstance(la2, float):
                    Latype2, a2 = "Square", la2
                else:
                    Latype2, a2 = Latype1, a1

                if isinstance(d_intra2, list) and len(d_intra2) == self.natom - 1:
                    d_a2 = d_intra2
                else:
                    d_a2 = d_a1
        else:
            elem2, Latype2, a2, d_a2 = elem1, Latype1, a1, d_a1
        return elem1, Latype1, a1, d_a1, elem2, Latype2, a2, d_a2

    def set_formula_w(self):
        if self.elements[0] != self.elements[2] and self.elements[3] != self.elements[5]:
            self.seqtype = 4
            self.formula_w = self.elements[1] + self.elements[0] + self.elements[2] + self.elements[4] + self.elements[
                3] + self.elements[-1]
        elif self.elements[0] == self.elements[2] and self.elements[3] == self.elements[5]:
            self.seqtype = 1
            self.formula_w = self.elements[1] + self.elements[0] + str(2) + self.elements[4] + self.elements[3] + str(2)
        elif self.elements[0] != self.elements[2] and self.elements[3] == self.elements[5]:
            self.seqtype = 2
            self.formula_w = self.elements[1] + self.elements[0] + self.elements[2] + self.elements[4] + self.elements[
                3] + str(2)
        else:
            self.seqtype = 3
            self.formula_w = self.elements[1] + self.elements[0] + str(2) + self.elements[4] + self.elements[3] + \
                             self.elements[5]
        if self.savenamemode == 1:
            self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + str(self.natom)

    def ini_positions(self):
        in_pos_half = 0.5 * self.la
        self.inplane_pos = {
            "AA": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [0.0, in_pos_half], [0.0, 0.0],
                   [in_pos_half, 0.0]],
            "AA2": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [in_pos_half, 0.0], [0.0, 0.0],
                    [0.0, in_pos_half]],
            "AB": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [0.0, in_pos_half], [in_pos_half, in_pos_half],
                   [in_pos_half, 0.0]],
            "AB2": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [in_pos_half, 0.0], [in_pos_half, in_pos_half],
                    [0.0, in_pos_half]],
            "AC": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [0.0, 0.0], [in_pos_half, 0.0],
                   [in_pos_half, in_pos_half]],
            "AC2": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [in_pos_half, in_pos_half], [in_pos_half, 0.0],
                    [0.0, 0.0]],
            "AC3": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [in_pos_half, in_pos_half], [0.0, in_pos_half],
                    [0.0, 0.0]],
            "AC4": [[0.0, in_pos_half], [0.0, 0.0], [in_pos_half, 0.0], [0.0, 0.0], [0.0, in_pos_half],
                    [in_pos_half, in_pos_half]],
        }

        self.outplane_pos = [self.lv / 2 + self.d_inter / 2 + sum(self.d_intra1),
                             self.lv / 2 + self.d_inter / 2 + self.d_intra1[-1],
                             self.lv / 2 + self.d_inter / 2,
                             self.lv / 2 - self.d_inter / 2,
                             self.lv / 2 - self.d_inter / 2 - self.d_intra2[0],
                             self.lv / 2 - self.d_inter / 2 - sum(self.d_intra2)]

        def pos_gen(stacktype):
            pos_ = []
            for i in range(6):
                pos_.append(np.append(self.inplane_pos[stacktype][i], self.outplane_pos[i]))
            return pos_

        self.stackpos = {}
        for postype_ in self.inplane_pos.keys():
            self.stackpos[postype_] = pos_gen(postype_)
        if self.hetero:
            self.stackmodes = list(self.stackpos.keys())
        else:
            self.stackmodes = list(self.stackpos.keys())[:-1]

        self.stackelements = {}
        if self.seqtype == 1:
            elem0 = self.elements[:]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0}
        elif self.seqtype == 2:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[2] = elem1[2], elem1[0]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1}
            self.stackelements["AA"] = {"cord1": elem0}
        elif self.seqtype == 3:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[3], elem1[5] = elem1[5], elem1[3]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1}
            self.stackelements["AA"] = {"cord1": elem0}
        else:
            elem0 = self.elements[:]
            elem1 = elem0[:]
            elem1[0], elem1[2] = elem1[2], elem1[0]
            elem2 = elem0[:]
            elem2[3], elem2[5] = elem2[5], elem2[3]
            elem3 = elem0[:]
            elem3[0], elem3[2], elem3[3], elem3[5] = elem3[2], elem3[0], elem3[5], elem3[3]
            for stack_mode in self.stackmodes:
                self.stackelements[stack_mode] = {"cord1": elem0, "cord2": elem1, "cord3": elem2, "cord4": elem3}

    def WritePOSCAR(self):
        """
        Write the POSCAR file for the given stacking pattern.
        """
        for stackmode in self.stackmodes:
            for cord_ in self.stackelements[stackmode]:
                element_pos_dict = {}
                pos_ = self.stackpos[stackmode]
                for i, e in enumerate(self.stackelements[stackmode][cord_]):
                    if e not in element_pos_dict:
                        element_pos_dict[e] = {}
                        element_pos_dict[e]['count'] = 0
                        element_pos_dict[e]['pos'] = []
                    element_pos_dict[e]['count'] += 1
                    element_pos_dict[e]['pos'].append(pos_[i])

                element_pos_dict = dict(sorted(element_pos_dict.items(), key=lambda x: self.elements.index(x[0])))

                if self.savenamemode == 1:
                    if not os.path.exists(os.path.join(self.savepath, self.formula_w, stackmode, cord_)):
                        os.makedirs(os.path.join(self.savepath, self.formula_w, stackmode, cord_))
                    savefile = os.path.join(self.savepath, self.formula_w, stackmode, cord_, f"POSCAR")
                else:
                    if not os.path.exists(self.savepath):
                        os.makedirs(self.savepath)
                    savefile = os.path.join(self.savepath,
                                            f"{self.formula_w}_{stackmode}_{cord_}-{self.mismatch:.2f}%-{str(self.natom)}-POSCAR")

                with open(savefile, 'w') as f:
                    f.write(f"{self.formula_w}-{stackmode}-{cord_}\n")
                    f.write("1.0\n")
                    f.write(f"{self.lattice_w[0][0]:.6f} {self.lattice_w[0][1]:.6f} {self.lattice_w[0][2]:.6f}\n")
                    f.write(f"{self.lattice_w[1][0]:.6f} {self.lattice_w[1][1]:.6f} {self.lattice_w[1][2]:.6f}\n")
                    f.write(f"{self.lattice_w[2][0]:.6f} {self.lattice_w[2][1]:.6f} {self.lattice_w[2][2]:.6f}\n")
                    for elem_ in element_pos_dict.keys():
                        f.write(f"{elem_} ")
                    f.write("\n")
                    for elem_ in element_pos_dict.keys():
                        f.write(f"{element_pos_dict[elem_]['count']} ")
                    f.write("\n")
                    f.write("Cartesian\n")
                    for elem_ in element_pos_dict.keys():
                        for pos_ in element_pos_dict[elem_]['pos']:
                            f.write(f"{pos_[0]:.6f} {pos_[1]:.6f} {pos_[2]:.6f}\n")

# class MNXYBilayerSquare(TMDHBilayerSquare):
#     def ini_atoms(self):
#         self.natom = 4
#
#     def set_formula_w(self):
#         if (self.elements[0] != self.elements[3] or self.elements[1] != self.elements[2]) and \
#                 (self.elements[4] != self.elements[-1] or self.elements[5] != self.elements[-2]):
#             self.seqtype = 4
#             self.formula_w = self.elements[1] + self.elements[2] + self.elements[0] + self.elements[3] + \
#                              self.elements[5] + self.elements[6] + self.elements[4] + self.elements[-1]
#         elif self.elements[0] == self.elements[3] and self.elements[1] == self.elements[2] and self.elements[4] == \
#                 self.elements[-1] and self.elements[5] == self.elements[-2]:
#             self.seqtype = 1
#             self.formula_w = self.elements[1] + str(2) + self.elements[0] + str(2) + self.elements[5] + str(2) \
#                              + self.elements[4] + str(2)
#         elif (self.elements[0] != self.elements[3] or self.elements[1] != self.elements[2]) and self.elements[4] == \
#                 self.elements[-1] and self.elements[5] == self.elements[-2]:
#             self.seqtype = 2
#             self.formula_w = self.elements[1] + self.elements[2] + self.elements[0] + self.elements[3] \
#                              + self.elements[5] + str(2) + self.elements[4] + str(2)
#         else:
#             self.seqtype = 3
#             self.formula_w = self.elements[1] + str(2) + self.elements[0] + str(2) + self.elements[5] + \
#                              self.elements[6] + self.elements[4] + self.elements[-1]
#         if self.savenamemode == 1:
#             self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + str(self.natom)

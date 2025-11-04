from pymatgen.core import Structure
import os
import warnings
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from collections import Counter


# def revPOS(st):
#     lattice = st.lattice.matrix
#     elements = st.species
#     positions = st.frac_coords
#     positions[:, 2] = -positions[:, 2]
#     rev_st = Structure(lattice=lattice,
#                        species=elements,
#                        coords=positions,
#                        coords_are_cartesian=False)
#     return rev_st


class BilayerAuto:
    def __init__(self, st1, st2=None, d_inter=3.3, lv=30.0, savepath=None, savenamemode=1, mismatch_threshold=5.0,
                 hetero_set=False,
                 stackmode="AA", overwrite=True, symprec=0.01, angle_tolerance=5, similar_threshold=0.01):
        """
        Constructing 2D bilayer structures.

        st1: list or str
            The atoms in the first layer or the POSCAR file of the first layer.
            If st1 is a list, note that the sequence of atoms in the list is ordered by height in Z direction.
            If st1 is a str, it is the path to the POSCAR file of the first layer, and a1, d_intra1 will be ignored.
        st2: list or str
            The atoms in the second layer or the POSCAR file of the second layer.
            If not provided, homogeneous bilayer will be constructed.
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
        symprec: float
            The symmetry precision for space group analysis.
        angle_tolerance: float
            The angle tolerance for space group analysis.
        similar_threshold: float
            The threshold to determine whether two structures are similar.
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
        self.overwrite = overwrite
        self.similar_threshold = similar_threshold
        self.matcher = StructureMatcher()
        self.hetero = False

        assert isinstance(st1, str), "The first layer st1 should be a POSCAR file path"
        self.mono1 = Structure.from_file(st1)
        st1_spg_analyzer = SpacegroupAnalyzer(self.mono1, symprec=symprec, angle_tolerance=angle_tolerance)
        self.Latype = st1_spg_analyzer.get_lattice_type()

        if st2 is not None:
            assert isinstance(st2, str), "The second layer st2 should be a POSCAR file path"
            self.mono2 = Structure.from_file(st2)
            st2_spg_analyzer = SpacegroupAnalyzer(self.mono2, symprec=symprec, angle_tolerance=angle_tolerance)
            st2_lattice_type = st2_spg_analyzer.get_lattice_type()
            assert self.Latype == st2_lattice_type, "The lattice types of the two layers are different."
            are_similar = self.matcher.fit(self.mono1, self.mono2)
            if are_similar:
                score = self.matcher.get_rms_dist(self.mono1, self.mono2)
                if score[1] < self.similar_threshold:
                    self.hetero = False
                else:
                    self.hetero = True
            else:
                self.hetero = True

        self.lv = lv
        self.stackmode = stackmode
        if hetero_set:
            self.hetero = True

        self.pos1_xy = np.array([site.frac_coords[:2] for site in self.mono1.sites])
        self.pos1_z = np.array([site.coords[2] for site in self.mono1.sites])

        if not self.hetero:
            elem = [str(site.specie) for site in self.mono1.sites]
            elem2 = elem
            self.lattice_w = self.mono1.lattice.matrix.copy()
            self.lattice_w[2, 2] = self.lv
            self.pos2_xy = self.pos1_xy.copy()
            self.pos2_z = self.pos1_z.copy()
        else:
            elem = [str(site.specie) for site in self.mono1.sites]
            elem2 = [str(site.specie) for site in self.mono2.sites]
            self.lattice_w = self.check_lattice_mismatch(self.mono1.lattice.a,
                                                         self.mono1.lattice.b,
                                                         self.mono1.lattice.angles[2],
                                                         self.mono2.lattice.a,
                                                         self.mono2.lattice.b)
            self.pos2_xy = np.array([site.frac_coords[:2] for site in self.mono2.sites])
            self.pos2_z = np.array([site.coords[2] for site in self.mono2.sites])

        self.elements = elem + elem2
        self.natom1 = len(elem)
        self.seq_cord = ["cord1"]

        revmono1_flag = False
        elem_reduce = self.reduce_element_pos(elem, self.pos1_z)
        if not self.hetero:
            if elem_reduce != elem_reduce[::-1]:
                self.seq_cord += ["cord2", "cord3"]
        else:
            if elem_reduce != elem_reduce[::-1]:
                revmono1_flag = True
            elem2_reduce = self.reduce_element_pos(elem2, self.pos2_z)
            if elem2_reduce != elem2_reduce[::-1]:
                if revmono1_flag:
                    self.seq_cord += ["cord2", "cord3", "cord4"]
                else:
                    self.seq_cord += ["cord3"]
            else:
                if revmono1_flag:
                    self.seq_cord += ["cord2"]

        count = Counter(self.elements)
        self.formula_w = ''.join(f"{element}{count[element] if count[element] > 1 else ''}" for element in count)
        if self.savenamemode == 1:
            self.formula_w += "_" + f"{self.mismatch:.2f}%" + "_" + f"{self.natom1}"

        self.d_inter = d_inter
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

    def check_lattice_mismatch(self, la1, lb1, angle1, la2, lb2):
        mislattice_a = 2 * abs(la1 - la2) / (la1 + la2) * 100
        mislattice_b = 2 * abs(lb1 - lb2) / (lb1 + lb2) * 100
        self.mismatch = max(mislattice_a, mislattice_b)
        la_mean = (la1 + la2) / 2
        lb_mean = (lb1 + lb2) / 2
        if angle1 > 90:
            angle1 = 180 - angle1
            angle1 = np.cos(np.radians(angle1))
            return np.array(
                [[la_mean, 0, 0], [- lb_mean * angle1, lb_mean * np.sqrt(1 - angle1 ** 2), 0],
                 [0, 0, self.lv]])
        else:
            angle1 = np.cos(np.radians(angle1))
            return np.array(
                [[la_mean, 0, 0], [lb_mean * angle1, lb_mean * np.sqrt(1 - angle1 ** 2), 0], [0, 0, self.lv]])

    def reduce_element_pos(self, elements, pos_z):
        """
        Reduces the elements and positions to the elements that are not in the same plane.
        After reduction, some stacking orders may be ignored.
        """
        elem = []
        z_coords_diff1 = np.diff(pos_z)
        for i, z_diff in enumerate(z_coords_diff1):
            if z_diff > 0.01:
                if len(elem) == 0:
                    elem.append(elements[i])
                    elem.append(elements[i + 1])
                else:
                    elem.append(elements[i + 1])
        return elem

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
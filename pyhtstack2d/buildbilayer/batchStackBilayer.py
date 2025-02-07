import os

from .stackBilayer import Monolayer, Bilayer, TMDHBilayer, N2Bilayer, MNXYBilayer, \
    N2TMDHbilayer, N2MNXYBilayer, TMDHBilayerSquare


gen_modes = {"bilayer": Bilayer, "tmdh": TMDHBilayer, "mnxy": MNXYBilayer, "n2": N2Bilayer,
             "n2tmdh": N2TMDHbilayer, "n2mnxy": N2MNXYBilayer, "tmdhsquare": TMDHBilayerSquare}
modes_kwargs = {"bilayer": Bilayer.__init__.__code__.co_varnames[3:Bilayer.__init__.__code__.co_argcount],
                "tmdh": TMDHBilayer.__init__.__code__.co_varnames[3:TMDHBilayer.__init__.__code__.co_argcount],
                "mnxy": MNXYBilayer.__init__.__code__.co_varnames[3:MNXYBilayer.__init__.__code__.co_argcount],
                "n2": N2Bilayer.__init__.__code__.co_varnames[3:N2Bilayer.__init__.__code__.co_argcount],
                "n2tmdh": N2TMDHbilayer.__init__.__code__.co_varnames[3:N2TMDHbilayer.__init__.__code__.co_argcount],
                "n2mnxy": N2MNXYBilayer.__init__.__code__.co_varnames[3:N2MNXYBilayer.__init__.__code__.co_argcount],
                "tmdhsquare": TMDHBilayerSquare.__init__.__code__.co_varnames[3:TMDHBilayerSquare.__init__.__code__.co_argcount]
                }


class GenBiLayer:
    def __init__(self, pos_dir, pos_dir2=None, la_mismatch=5.0, lb_mismatch=None, homo=False, genmode="bilayer",
                 **kwargs):
        """
        Batch stacking of bilayer homogeneous/heterogeneous structures.

        pos_dir: str or list
            The directory containing the POSCAR files of the structures to be stacked.
            If pos_dir is a string, all the POSCAR files in the directory will be stacked.
            If pos_dir is a list, only the POSCAR files in the list will be stacked.
        pos_dir2: str or list
            If pos_dir2 is not None, iterate through all combinations of pos_dir and pos_dir2.
            Note that the internal combinations of pos_dir (pos_dir2) are not considered in this case, and parameters homo is invalid.
        la_mismatch: float/int or list/tuple
            The maximum lattice mismatch allowed between the two layers, the default value is 5.0.
            If la_mismatch is float or int, stacking will not be performed if the lattice mismatch is greater than this value.
            If la_mismatch is list or tuple, the lattice mismatch should be within la_mismatch.
            For example, la_mismatch=[3.0, 5.0], which means that only the lattice mismatch is in [3.0, 5.0], stacking will be performed.
        lb_mismatch: float
            When the in-plane lattice vectors are different, the second lattice mismatch should be specified.
            The maximum allowed lattice mismatch between the two layers.
            If lb_mismatch is float or int, stacking will not be performed if the lattice mismatch is greater than this value.
            If lb_mismatch is list or tuple, the lattice mismatch should be within lb_mismatch.
        homo: bool
            If True, both homogeneous and heterogeneous structures will be stacked.
            If False, only heterogeneous structures will be stacked.
        genmode: str
            The mode of stacking. The available modes are:
            "bilayer": Stacking of two bilayers.
            "tmdbilayer": Stacking of a bilayer on a TMD monolayer.
            "n2bilayer": Stacking of a bilayer on a N2 monolayer.
        **kwargs: pos_dict
            The keyword arguments for the stacking mode.
            Use the GenBiLayer.get_kwargs() method to get the available keyword arguments for the specified mode.
        """
        assert genmode.lower() in gen_modes.keys(), f"genmode should be one of {gen_modes.keys()}"

        self.single_dir = True
        self.pos_obj, self.la, self.lb = self.get_pos_inf(pos_dir)
        if pos_dir2:
            self.pos_obj2, self.la2, self.lb2 = self.get_pos_inf(pos_dir2)
            self.single_dir = False

        self.genmode = genmode.lower()
        self.la_mismatch = la_mismatch
        self.lb_mismatch = lb_mismatch
        if isinstance(self.la_mismatch, (int, float)):
            self.la_upper = self.la_mismatch
            self.la_lower = 0.0
        elif isinstance(self.la_mismatch, (list, tuple)):
            assert len(self.la_mismatch) == 2, "la_mismatch should be a float or a list/tuple of two floats."
            self.la_upper = self.la_mismatch[1] if self.la_mismatch[1] > self.la_mismatch[0] else self.la_mismatch[0]
            self.la_lower = self.la_mismatch[0] if self.la_mismatch[1] > self.la_mismatch[0] else self.la_mismatch[1]
        else:
            raise ValueError("la_mismatch should be a float or a list/tuple of two floats.")
        if self.lb_mismatch:
            if isinstance(self.lb_mismatch, (int, float)):
                self.lb_upper = self.lb_mismatch
                self.lb_lower = 0.0
            elif isinstance(self.lb_mismatch, (list, tuple)):
                assert len(self.lb_mismatch) == 2, "lb_mismatch should be a float or a list/tuple of two floats."
                self.lb_upper = self.lb_mismatch[1] if self.lb_mismatch[1] > self.lb_mismatch[0] else self.lb_mismatch[0]
                self.lb_lower = self.lb_mismatch[0] if self.lb_mismatch[1] > self.lb_mismatch[0] else self.lb_mismatch[1]
            else:
                raise ValueError("lb_mismatch should be a float or a list/tuple of two floats.")
        self.homo = homo

        for kwarg in kwargs.keys():
            assert kwarg in modes_kwargs[self.genmode], \
                f"kwarg of \"{self.genmode}\" mode should be one of {modes_kwargs[self.genmode]}"
        self.kwargs = kwargs
        if self.single_dir:
            self.match_pos_dict = self.match_pos()
        else:
            self.match_pos_dict = self.match_pos2()

    def get_pos_inf(self, pos_dir):
        pos_obj = []
        la = []
        lb = []
        if isinstance(pos_dir, str):
            for pos_file in os.listdir(pos_dir):
                momolayer = Monolayer(os.path.join(pos_dir, pos_file))
                pos_obj.append(momolayer)
                la.append(momolayer.a)
                lb.append(momolayer.b)
        elif isinstance(pos_dir, list):
            for pos_file in pos_dir:
                momolayer = Monolayer(pos_file)
                pos_obj.append(momolayer)
                la.append(momolayer.a)
                lb.append(momolayer.b)
        else:
            raise ValueError("pos_obj should be a string or a list of strings")
        return pos_obj, la, lb

    def match_pos(self):
        match_pos_dicr = {}
        if self.homo:
            index_add = 0
        else:
            index_add = 1

        for i_index, la_i in enumerate(self.la):
            match_pos_dicr[i_index] = []
            lb_i = self.lb[i_index]
            for j_index, la_j in enumerate(self.la[i_index + index_add:]):
                lb_j = self.lb[i_index + index_add + j_index]
                delta_la = 2 * abs(la_i - la_j) / (la_i + la_j) * 100
                delta_lb = 2 * abs(lb_i - lb_j) / (lb_i + lb_j) * 100
                if self.lb_mismatch is None and self.la_lower <= delta_la <= self.la_upper:
                    match_pos_dicr[i_index].append(i_index + index_add + j_index)
                elif self.lb_mismatch and self.la_lower <= delta_la <= self.la_upper and self.lb_lower <= delta_lb <= self.lb_upper:
                    match_pos_dicr[i_index].append(i_index + index_add + j_index)
                else:
                    continue
        return match_pos_dicr

    def match_pos2(self):
        match_pos_dicr = {}
        for i_index, la_i in enumerate(self.la):
            match_pos_dicr[i_index] = []
            lb_i = self.lb[i_index]
            for j_index, la_j in enumerate(self.la2):
                lb_j = self.lb2[j_index]
                delta_la = 2 * abs(la_i - la_j) / (la_i + la_j) * 100
                delta_lb = 2 * abs(lb_i - lb_j) / (lb_i + lb_j) * 100
                if self.lb_mismatch is None and self.la_lower <= delta_la <= self.la_upper:
                    match_pos_dicr[i_index].append(j_index)
                elif self.lb_mismatch and self.la_lower <= delta_la <= self.la_upper and self.lb_lower <= delta_lb <= self.lb_upper:
                    match_pos_dicr[i_index].append(j_index)
                else:
                    continue
        return match_pos_dicr

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
                    gen_modes[self.genmode](st1=mono1, st2=mono2, lv=lv, **self.kwargs).WritePOSCAR()
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
                    gen_modes[self.genmode](st1=mono1, st2=mono2, lv=lv, **self.kwargs).WritePOSCAR()

        savedir = gen_modes[self.genmode](self.pos_obj[0], **self.kwargs).savepath
        print(f"Stacking of bilayers is completed. The POSCAR files are saved in the {savedir} directory.")

    @staticmethod
    def get_kwargs():
        for mode_ in modes_kwargs.keys():
            print(f"Specified kwargs for the \"{mode_}\" mode: \n {modes_kwargs[mode_]}\n")

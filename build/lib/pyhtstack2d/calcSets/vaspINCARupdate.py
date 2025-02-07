import os

import numpy as np
from .Pu import INCARPu


class UpdateINCAR:
    def __init__(self, materiallist="material.txt", dirprefix=None, incarbasic=None, pU=True, ismag=True,
                 collateincar=False, **kwargs):
        """
        Update INCAR file with new parameters or from ICNAR-basic file.

        materiallist: list of materials
        dirprefix: prefix of the directories containing the materials
        incarbasic: path to the basic INCAR file or dictionary of parameters
        pU: if True, add U setting to the INCAR file
        collateincar: if True, collate INCAR content after updating. Sometimes the generated INCAR will have duplicate parameter lines, this parameter can be used to rearrange the contents of the INCAR
        ismag: whether the calculation is magnetic, if True, it will search the directories FM, AFM1, AFM2 etc.

        """
        assert os.path.exists(materiallist), "Material list not found!"
        assert incarbasic is not None, "INCAR-basic file or dictionary of parameters not provided!"

        with open(materiallist, 'r') as file:
            self.materiallist = [line.strip() for line in file]
        self.dirprefix = dirprefix or ""
        self.incarbasic = incarbasic
        self.pU = pU
        self.ismag = ismag
        self.collateincar = collateincar

        self.kwargs = kwargs
        self.puset = {"u_setting": None, "magatomlist": None, "openmixparam": True, "mixparam": None}
        for key in self.puset.keys():
            if key in self.kwargs:
                self.puset[key] = self.kwargs

    def parse_icnar(self, incarbasic):
        """Parse ICNAR file to extract parameters."""
        params = {}
        with open(incarbasic, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                key, value = line.split('=')
                params[key.strip()] = value.strip()
        return params

    def update_incar(self, incarpath, newparams):
        """Update INCAR file with new parameters."""
        incar_lines = []
        with open(incarpath, 'r') as file:
            for line in file:
                updated = False
                for key, value in newparams.items():
                    if line.startswith(key):
                        incar_lines.append(f"{key} = {value}\n")
                        updated = True
                        break
                if not updated:
                    incar_lines.append(line)

        with open(incarpath, 'w') as file:
            for line in incar_lines:
                file.write(line)
            for key, value in newparams.items():
                if f"{key} = " not in "".join(incar_lines):
                    file.write(f"{key} = {value}\n")

    def collate_incar(self, incar_file):
        with open(incar_file, 'r') as file:
            lines = file.readlines()
        parameters = {}
        for line in lines:
            stripped_line = line.strip()
            if stripped_line and not stripped_line.startswith('#'):
                if '=' in stripped_line:
                    key, value = stripped_line.split('=', 1)
                    parameters[key.strip()] = value.strip()
        with open(incar_file, 'w') as file:
            for key, value in parameters.items():
                file.write(f"{key} = {value}\n")

    def updateIncarfiles(self):
        """Update all INCAR files in directories with parameters from ICNAR file."""
        if isinstance(self.incarbasic, str):
            assert os.path.exists(self.incarbasic), "INCAR-basic file not found!"
            incarbasicpara = self.parse_icnar(self.incarbasic)
        else:
            incarbasicpara = self.incarbasic

        for incarpath in self.materiallist:
            incar_file = os.path.join(incarpath, self.dirprefix, 'INCAR')
            if os.path.exists(incar_file) and not os.path.islink(incar_file):
                newparams = incarbasicpara.copy()
                if self.pU:
                    poscar_file = os.path.join(incarpath, 'POSCAR')
                    if os.path.exists(poscar_file):
                        elements = np.loadtxt(poscar_file, skiprows=5, max_rows=1, dtype=str)
                        upara = INCARPu(elements, incar=None, u_setting=self.puset["u_setting"], magatomlist=self.puset["magatomlist"],
                                        openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                        for key, value in upara.items():
                            if isinstance(value, list):
                                upara[key] = " ".join(map(str, value))
                        newparams.update(upara)
                self.update_incar(incar_file, newparams)
                if self.collateincar:
                    self.collate_incar(incar_file)
                # print(f"Updated INCAR in {incarpath} with new parameters.")

            if self.ismag:
                for dir_ in ["FM", "AFM1", "AFM2", "AFM3"]:
                    incar_file = os.path.join(incarpath, dir_, 'INCAR')
                    if os.path.exists(incar_file) and not os.path.islink(incar_file):
                        newparams = incarbasicpara.copy()
                        if self.pU:
                            poscar_file = os.path.join(incarpath, dir_, 'POSCAR')
                            if os.path.exists(poscar_file):
                                elements = np.loadtxt(poscar_file, skiprows=5, max_rows=1, dtype=str)
                                upara = INCARPu(elements, incar=None, u_setting=self.puset["u_setting"], magatomlist=self.puset["magatomlist"],
                                        openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                                for key, value in upara.items():
                                    if isinstance(value, list):
                                        upara[key] = " ".join(map(str, value))
                                newparams.update(upara)
                        self.update_incar(incar_file, newparams)
                        if self.collateincar:
                            self.collate_incar(incar_file)
                        # print(f"Updated INCAR in {incarpath}/{dir_} with new parameters.")

    def genIncarfiles(self):
        """Generate new INCAR files in directories from INCAR-basic file."""
        if isinstance(self.incarbasic, dict):
            with open("INCAR-basic", 'w') as file:
                for key, value in self.incarbasic.items():
                    file.write(f"{key} = {value}\n")
            incarbasicfile = "INCAR-basic"
        else:
            incarbasicfile = self.incarbasic

        for incarpath in self.materiallist:
            incar_file = os.path.join(incarpath, self.dirprefix, 'INCAR')
            if not os.path.islink(incar_file) or self.pU:
                os.system(f"cp {incarbasicfile} {incar_file}")
                if self.pU:
                    poscar_file = os.path.join(incarpath, 'POSCAR')
                    if os.path.exists(poscar_file):
                        elements = np.loadtxt(poscar_file, skiprows=5, max_rows=1, dtype=str)
                        INCARPu(elements, incar=incar_file, u_setting=self.puset["u_setting"], magatomlist=self.puset["magatomlist"],
                                        openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                if self.collateincar:
                    self.collate_incar(incar_file)
                # print(f"Generated INCAR in {incarpath}.")

            if self.ismag:
                for dir_ in ["FM", "AFM1", "AFM2", "AFM3"]:
                    incar_file = os.path.join(incarpath, dir_, 'INCAR')
                    if os.path.exists(os.path.join(incarpath, dir_)):
                        os.system(f"cp {incarbasicfile} {incar_file}")
                        if self.pU:
                            poscar_file = os.path.join(incarpath, dir_, 'POSCAR')
                            if os.path.exists(poscar_file):
                                elements = np.loadtxt(poscar_file, skiprows=5, max_rows=1, dtype=str)
                                INCARPu(elements, incar=incar_file, u_setting=self.puset["u_setting"], magatomlist=self.puset["magatomlist"],
                                        openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                        if self.collateincar:
                            self.collate_incar(incar_file)
                        # print(f"Generated INCAR in {incarpath}/{dir_}.")

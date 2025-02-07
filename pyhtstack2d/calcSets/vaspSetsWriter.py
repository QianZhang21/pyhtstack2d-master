"""
This module is used to set up the input files for VASP calculations, in particular some basic parameter settings
for the ground state calculations (relax, scf, band, etc.) for 2D materials.

According to the common settings configuration, users can modify the json file themselves.
INCARPara.json: The common settings for INCAR file.
pmap.json: The POTCAR map for the elements.
maginit.json: The common settings for the initial magnetic moment.
"""
import os
import json
import warnings

import numpy as np
# import shutil
from shutil import copyfile

from .IncarVePu import INCARSetVePu
from .Pu import INCARPu

INCARParaFile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'INCARPara.json')
INCARParaFile_reset = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'INCARPara_reset.json')
with open(INCARParaFile, "r") as f:
    INCAR_para = json.load(f)

POT_PATH = os.path.join(os.path.expanduser("~"), ".config", ".PyHTStack2D.json")

pmappath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pmap.json')
with open(pmappath, "r") as f:
    pmap = json.load(f)

maginitpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'maginit.json')
with open(maginitpath, "r") as f:
    maginit = json.load(f)


class IncParaSet:
    def __init__(self, parameter=None):
        """
        Update parameter settings commonly used by users in the INCARPara.json file.

        parameter: pos_dict
            The parameter settings to be updated.
        """
        self.parameter = parameter
        self.INCAR_para = INCAR_para

    def update(self):
        """
        Update the parameter settings in the INCARPara.json file.
        """
        if self.parameter:
            self.INCAR_para.update(self.parameter)
            print("Update the INCARPara.json file.")
            with open(INCARParaFile, "w") as f:
                json.dump(self.INCAR_para, f, indent=4)

    def get(self):
        """
        Get the specific parameter settings in the INCARPara.json file.
        """
        return self.INCAR_para

    def get_incar_set(self):
        """
        Get the keys of the parameter settings in the INCARPara.json file.
        """
        opts = list(INCAR_para.keys())
        opts.remove("LISTPara")
        print(f"The available calculation modes are as follows: {opts}")
        return opts

    def reset(self):
        """
        Reset the parameter settings in the INCARPara.json file.
        """
        if os.path.exists(INCARParaFile_reset):
            copyfile(INCARParaFile_reset, INCARParaFile)
            print("Reset the INCARPara.json file.")


class StructureFromFile:
    def __init__(self, stfile):
        """
        Set the structure file for VASP calculation.

        stfile: str
            Path to the POSCAR file.
        """
        self.StFile = stfile
        self.uid = os.path.basename(self.StFile).split('-')[0]
        self.latticeA = None
        self.latticeB = None
        self.latticeC = None
        self.lattice = None
        self.elements = None
        self.elemnum = None
        self.poscarType = None
        self.positions = None
        self.getInfo()

    def getInfo(self):
        """
        Get the information of the structure file.
        """
        if not os.path.exists(self.StFile):
            raise FileNotFoundError("The POSCAR file does not exist.")

        with open(self.StFile, "r") as f:
            lines = f.readlines()

        self.lattice = np.array([list(map(float, lines[i].split())) for i in range(2, 5)])
        self.latticeA = self.lattice[0]
        self.latticeB = self.lattice[1]
        self.latticeC = self.lattice[2]
        self.elements = list(map(str, lines[5].split()))
        self.elemnum = list(map(int, lines[6].split()))
        num_atoms = sum(self.elemnum)
        self.poscarType = lines[7].strip().lower()
        self.positions = np.array([list(map(float, line.split()[:3])) for line in lines[8:8 + num_atoms]])


class INCARSet:
    def __init__(self, savepath, stfile=None, inputOpt="scf", user_incar_settings=None, ini_magmom=False,
                 magmom_setting=None, nonmagmom=0):
        """
        Set INCAR parameters for VASP calculation.

        savepath: str
            Path to save the INCAR file.
        stfile: str, default=None
            Path to the POSCAR file.
        inputOpt: str, default="scf"
            The calculation mode of VASP. The basic modes available are 'relax', 'scf', 'band', 'optic', 'eps' (electrostatic potential).
            Some additional suffixes are available: '-hse06' (HSE06 hybrid functional), '-pu' (DFT + U),
            '-dip' (Dipole correction), <'-d2', '-d3', '-d3bj', '-optpbe', '-optb88', '-optb86b'> (vdW correction).
            For example, 'scf-d3' means the DFT-D3 no-damping correction static calculation.
            Through the get_opt() method, the user can obtain the available calculation modes, e.g., INCARSet("").get_opt().
        user_incar_settings: pos_dict
            User-defined INCAR settings.
            Generally, 99% of tasks require specialization of the INCAR file by user_incar_settings.
            For example, user_incar_settings = {"NPAR": 4, "EDIFF": 1e-5, , "IVDW": 12,
            "LDAUU": [3.9, 0], "LDAUJ": [1.1, 0]}
            Note that if the format of "MAGMOM" in user_incar_settings is NIONS*1.0, the elements of the list should
            be str, such as user_incar_settings={"MAGMOM": ["2*5.0", "18*0.6"]}.
            "ISPIN" is set to 2 if usser_incar_settings contains "MAGMOM".
        ini_magmom: bool, default=False
            Whether to set the initial magnetic moment for the calculation.
        magmom_setting: dict, default=None
            The initial magnetic moment settings for the calculation.
        nonmagmom: float, default=0
            The initial magnetic moment for non-magnetic atoms.
        """
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        self.savepath = os.path.join(savepath, "INCAR")
        self.stfile = stfile
        self.cal_mode = inputOpt.split("-")
        self.base_incar = INCAR_para["scf"]
        # self.ispin = 1
        self.ini_magmom = ini_magmom
        self.magmom_setting = magmom_setting
        self.nonmagmom = nonmagmom
        self.user_incar_settings = user_incar_settings
        self.check_mode_para()
        self.incar_write = {}
        self.elements = None
        self.elemnum = None

    def check_mode_para(self):
        """
        Check the INCAR modes and the List parameters in user_incar_settings.
        """
        for mode in self.cal_mode:
            if mode.lower() not in INCAR_para.keys() and mode.lower() != 'pu':
                av_k = list(INCAR_para.keys())
                av_k.remove("LISTPara")
                print("Only the following calculation modes are currently supported: ", av_k)
                raise KeyError(f"The calculation mode {mode} is not supported.")

        if self.user_incar_settings is not None:
            for key, value in self.user_incar_settings.items():
                if (isinstance(value, list) and key not in INCAR_para["LISTPara"]) or \
                        (not isinstance(value, list) and key in INCAR_para["LISTPara"]):
                    print("These parameters must be of type list: ", INCAR_para["LISTPara"])
                    raise KeyError(f"The parameter {key} may not be the correct type.")

    def incar_updata(self):
        """
        Update the INCAR parameters according to the calculation mode and user-defined settings.
        """
        self.incar_write.update(self.base_incar)
        for mode in self.cal_mode:
            if mode != 'pu':
                self.incar_write.update(INCAR_para[mode])
            else:
                if self.elements is None:
                    if isinstance(self.stfile, StructureFromFile):
                        pass
                    elif isinstance(self.stfile, str):
                        self.stfile = StructureFromFile(self.stfile)
                    elif self.stfile is None:
                        stf_ = os.path.dirname(self.savepath).replace('INCAR', 'POSACR')
                        self.stfile = StructureFromFile(stf_)
                    else:
                        raise ValueError("The stfile parameter is not correct.")
                    self.elements = self.stfile.elements
                pu = INCARPu(self.elements)
                self.incar_write.update(pu)

        if self.user_incar_settings is not None:
            self.incar_write.update(self.user_incar_settings)
            if "MAGMOM" in self.user_incar_settings.keys():
                self.incar_write["ISPIN"] = 2
                for mag_ in self.user_incar_settings["MAGMOM"]:
                    if isinstance(mag_, (float, int)):
                        warnings.warn("Please note the format of MAGMOM. "
                                      "If the format is NIONS*1.0, the elements of the list should be str, "
                                      "e.g. user_incar_settings={\"MAGMOM\": ['2*5.0', '18*0.6']}.")
                        break
        # self.ispin = self.incar_write["ISPIN"]

    def ini_magmom_set(self):
        """
        Set the initial magnetic moment for the calculation.
        """
        if self.elements is None or self.elemnum is None:
            if isinstance(self.stfile, StructureFromFile):
                pass
            elif isinstance(self.stfile, str):
                self.stfile = StructureFromFile(self.stfile)
            elif self.stfile is None:
                stf_ = os.path.dirname(self.savepath).replace('INCAR', 'POSACR')
                self.stfile = StructureFromFile(stf_)
            else:
                raise ValueError("The stfile parameter is not correct.")
            self.elements = self.stfile.elements
            self.elemnum = self.stfile.elemnum

        magmom = []
        for i, elem in enumerate(self.elements):
            if self.magmom_setting is not None and elem in self.magmom_setting.keys():
                magmom.extend([str(self.elemnum[i]) + '*' + str(self.magmom_setting[elem])])
            elif elem in maginit["MAGMOM"].keys():
                magmom.extend([str(self.elemnum[i]) + '*' + str(maginit["MAGMOM"][elem])])
            else:
                magmom.extend([str(self.elemnum[i]) + '*' + str(self.nonmagmom)])
                # if isinstance(self.ini_magmom_val, (float, int)):
                #     magmom.extend([str(self.elemnum[i]) + '*' + str(self.ini_magmom_val)])
                # elif isinstance(self.ini_magmom_val, dict):
                #     if elem in self.ini_magmom_val.keys():
                #         magmom.extend([str(self.elemnum[i]) + '*' + str(self.ini_magmom_val[elem])])
                #     else:  # The default value is 0 if the element is not in the pos_dict.
                #         magmom.extend([str(self.elemnum[i]) + '*' + str(0)])
        self.incar_write["MAGMOM"] = magmom
        self.incar_write["ISPIN"] = 2

    def writefile(self):
        """
        Write the INCAR file.
        """
        self.incar_updata()
        if self.ini_magmom and "MAGMOM" not in self.incar_write.keys():
            self.ini_magmom_set()

        list_para_w = {}
        for key, value in self.incar_write.items():
            if isinstance(value, list):
                list_para_w[key] = value
        for key in list_para_w.keys():
            self.incar_write.pop(key)

        with open(self.savepath, "w") as f:
            for key, value in self.incar_write.items():
                f.write(f"{key} = {value}\n")
            for key, value in list_para_w.items():
                f.write(f"{key} = ")
                for v in value:
                    f.write(f"{v} ")
                f.write("\n")

    def get_opt(self):
        """
        Get the optional INCAR parameters in INCARPara.json.
        """
        opts = list(INCAR_para.keys())
        opts.remove("LISTPara")
        print(f"The available calculation modes are as follows: {opts}")
        return opts


class MeshKPOINTSet:
    def __init__(self, savepath, kselfSet=None, StObj=None, kmeshrv=0, gamma=True, is2D=True):
        """
        Only support the regular mesh to select k-points.
        For band structure calculations, users can use the KPOINTS file generated by the VASPKIT package or set the k-points manually.

        kselfSet: list, default=None
            The customised list, such as [11, 11, 1].
        StObj: StructureFromFile object or str
            The object of the StructureFromFile class or the path of the POSCAR file.
        kmeshrv: float,
            The input Kmesh-resolved value (in units of 2*PI/Angstrom), as in the VASPKIT package.
            About the VASPKIT, see https://vaspkit.com
            Gamma-only: 0, Low accuracy: 0.06~0.04, Medium accuracy (general selection): 0.03~0.04,
            Fine accuracy: 0.02-0.01.
        gamma: bool, default=False
            Whether to use the gamma-centered k-point.
        """
        assert isinstance(kselfSet, (list, type(None))), "The kselfSet must be a list, whose length is 3."
        assert isinstance(StObj, (StructureFromFile, str)), "The StObj must be a StructureFromFile object " \
                                                            "or the path of POSCAR."
        assert isinstance(kmeshrv, (int, float)) and kmeshrv >= 0, "The kmeshrv must be a non-negative number."
        assert isinstance(gamma, bool), "The gamma must be a boolean."
        assert isinstance(is2D, bool), "The is2D must be a boolean."

        if not os.path.exists(savepath):
            os.makedirs(savepath)
        self.savepath = os.path.join(savepath, "KPOINTS")
        if isinstance(kselfSet, list) and len(kselfSet) == 3:
            self.KP = kselfSet
        else:
            if isinstance(StObj, StructureFromFile):
                self.StObj = StObj
            else:
                self.StObj = StructureFromFile(StObj)
            self.kmeshrv = kmeshrv
            self.gamma = gamma
            self.is2D = is2D
            self.volume = None
            self.Recip_a = self.StObj.latticeA
            self.Recip_b = self.StObj.latticeB
            self.Recip_c = self.StObj.latticeC
            self.poscarType = self.StObj.poscarType
            self.position = self.StObj.positions
            self.c_ = 0
            self.kmeshgen(is2D)

        self.ktype = "Gamma" if self.gamma else "Monkhorst-Pack"

    def kmeshgen(self, is2D=True):
        """
        Calculate the K-mesh according to the Kmesh-resolved value.
        """
        if self.kmeshrv == 0:
            self.KP = [1, 1, 1]
        else:
            a_ = self.StObj.latticeA
            b_ = self.StObj.latticeB
            c_ = self.StObj.latticeC

            self.volume = np.dot(a_, np.cross(b_, c_))
            if self.volume <= 0:
                raise ValueError("The volume of the cell is not positive.")

            self.Recip_a = np.linalg.norm(np.cross(b_, c_) / self.volume)
            self.Recip_b = np.linalg.norm(np.cross(c_, a_) / self.volume)
            self.Recip_c = np.linalg.norm(np.cross(a_, b_) / self.volume)

            if is2D:
                self.poscarType = self.StObj.poscarType
                self.position = self.StObj.positions
                self.c_ = np.linalg.norm(c_)

            dz = 0

            if self.is2D and self.position is not None:
                dz = np.max(self.position[:, 2]) - np.min(self.position[:, 2])
                if self.poscarType.lower().startswith('d'):
                    dz *= self.c_
                dz = abs(dz - self.c_)

            kx = round(self.Recip_a / self.kmeshrv)
            ky = round(self.Recip_b / self.kmeshrv)
            kz = round(self.Recip_c / self.kmeshrv) if not self.is2D else 1

            if kx < 1:
                kx = 1
            if ky < 1:
                ky = 1

            if self.is2D and dz < 10:
                warnings.warn("The K-mesh in the z-direction is set to 1. \n"
                              "The vacuum layer may be too thick. Please check the POSCAR file.")

            self.KP = [kx, ky, kz]

    def writefile(self):
        """
        Write the KPOINTS file.
        """
        with open(self.savepath, 'w') as f:
            f.write(f"K-Mesh, accuracy value is {self.kmeshrv}\n")
            f.write("0\n")
            f.write(self.ktype + "\n")
            f.write("  ".join(map(str, self.KP)))
            f.write("\n")
            f.write("0.0  0.0  0.0")


class BandKPOINTSet:
    def __init__(self, savepath, postype="H", angle=120):
        """
        The class for generating the KPOINTS file for the band structure calculation.
        Currently, only linear modes of k-path for hexagonal and rectangular lattices are supported for band calculations.

        savepath: str
            The path to save the KPOINTS file.
        postype: str, default="H"
            The type of the POSCAR file. The default value is "H" for hexagonal lattice.
            "R" for the rectangular lattice.
        angle: float, default=120
            The angle between the two lattice vectors for the hexagonal lattice.
        """
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        if postype.upper() not in ["H", "R"]:
            warnings.warn("The postype must be 'H' (hexagonal) or 'R' (rectangular). Now the postype is set to 'H'.")
            postype = "H"

        self.savepath = os.path.join(savepath, "KPOINTS")
        self.postype = postype.upper()
        self.angle = angle

    def writefile(self):
        """
        Write the KPOINTS file.
        """
        with open(self.savepath, 'w') as f:
            f.write("K-Path for band structure calculation\n")
            f.write("  60\nLine-Mode\nReciprocal\n")
            if self.postype == "H":
                if abs(self.angle-120) <= 1e-2 or abs(self.angle+0.5) <= 1e-2 :
                    f.write("  0.0000000000   0.0000000000   0.0000000000     GAMMA\n")
                    f.write("  0.3333333333   0.0000000000   0.0000000000     K\n\n")
                    f.write("  0.3333333333   0.0000000000   0.0000000000     K\n")
                    f.write("  0.5000000000   0.0000000000   0.0000000000     M\n\n")
                    f.write("  0.5000000000   0.0000000000   0.0000000000     M\n")
                    f.write("  0.0000000000   0.0000000000   0.0000000000     GAMMA\n")
                elif abs(self.angle-60) <= 1e-2 or abs(self.angle-0.5) <= 1e-2:
                    f.write("  0.0000000000   0.0000000000   0.0000000000     GAMMA\n")
                    f.write("  0.5000000000   0.5000000000   0.0000000000     M\n\n")
                    f.write("  0.5000000000   0.5000000000   0.0000000000     M\n")
                    f.write("  0.3333333333   0.6666666667   0.0000000000     K\n\n")
                    f.write("  0.3333333333   0.6666666667   0.0000000000     K\n")
                    f.write("  0.0000000000   0.0000000000   0.0000000000     GAMMA\n")
                else:
                    warnings.warn("The angle is not 120 or 60 degrees. The KPOINTS file is not generated.")
            else:
                f.write("  0.0   0.0   0.0     GAMMA\n")
                f.write("  0.5   0.0   0.0     X\n\n")
                f.write("  0.5   0.0   0.0     X\n")
                f.write("  0.5   0.5   0.0     S\n\n")
                f.write("  0.5   0.5   0.0     S\n")
                f.write("  0.0   0.5   0.0     Y\n\n")
                f.write("  0.0   0.5   0.0     Y\n")
                f.write("  0.0   0.0   0.0     GAMMA\n")


class POTCARSet:
    def __init__(self, savepath, elements, potpath=None):
        """
        Set POTCAR parameters for VASP calculation.

        elements: list or str
            The list of elements in the POSCAR file or the path of the POSCAR file.
        potpath: str, default=None
            The path of the POTCAR directory.
            For example, '/home/POTCAR_DIR/element/POTCAR', please set the potpath='/home/POTCAR_DIR'.
        """
        if isinstance(elements, str):
            assert os.path.exists(elements), "The POSCAR file does not exist."
            StObj = StructureFromFile(elements)
            elements = StObj.elements
        elif isinstance(elements, list):
            assert all(isinstance(el, str) for el in elements), "The elements must be a list of strings."
        else:
            raise ValueError("The elements must be a list of strings or the path of the POSCAR file.")
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        self.savepath = os.path.join(savepath, "POTCAR")
        self.elements = [pmap[el] for el in elements]
        self.potpath = potpath

        if self.potpath is None:
            assert os.path.exists(POT_PATH), "Please check if '~/.config/.PyHTStack2D.json' exists " \
                                             "and make sure the POTCAR path is added."
            potcar_config_ = json.load(open(POT_PATH, 'r'))
            self.potpath = potcar_config_["POTCAR_PATH"]
        else:
            self.potpath = potpath

    def writefile(self):
        """
        Get the POTCAR file according to the elements.
        """
        with open(self.savepath, 'w') as f:
            for element in self.elements:
                potcar_path = os.path.join(self.potpath, element, "POTCAR")
                if not os.path.exists(potcar_path):
                    raise FileNotFoundError(f"The POTCAR file of {element} does not exist.")
                with open(potcar_path, 'r') as potcar_r:
                    f.write(potcar_r.read())


class InputWriter:
    def __init__(self, stfile=None, posdir="POSCAR_dir", workdir="work", inputOpt='scf', overwrite=True,
                 user_incar_settings=None,
                 ini_magmom=False, magmom_setting=None, nonmagmom=0,
                 kselfSet=None, kmeshrv=0, gamma=True, postype=None, angle=120, is2D=True, potpath=None, is_print=False,
                 taskname=None):
        """
        Write the input files for VASP calculation.

        stfile: str or list
            The path of the structure file or the list of the structure file.
            When stfile is None, it will search the POSCAR files in the posdir.
        posdir: str, default='POSCAR_dir'
            The directory where the POSCAR files are located.
            For example, posdir = 'POSCAR_dir' means the POSCAR files are located in the 'POSCAR_dir' directory.
        workdir: str, default='work'
            The directory where the calculation files are saved.
        overwrite: bool, default=True
            Whether to overwrite the existing files.
        inputOpt: str, default='scf'
            The calculation mode, 'relax', 'scf', 'band', 'optic', 'eps' (electrostatic potential).
            Some additional suffixes are available: '-hse06' (HSE06 hybrid functional), '-pu' (DFT + U),
            '-dip' (Dipole correction), <'-d2', '-d3', '-d3bj', '-optpbe', '-optb88', '-optb86b'> (vdW correction).
            For example, 'scf-d3' means the DFT-D3 no-damping correction static calculation .
            If inputOpt='VePu', the INCAR file for the DFT+U and valance electrons calculation is generated.
        user_incar_settings: pos_dict, default=None
            The user-defined INCAR settings.
        ini_magmom: bool, default=False
            Whether to set the initial magnetic moment.
        magmom_setting: dict, default=None
            The magnetic moment settings.
        nonmagmom: float, default=0
            The non-magnetic moment.
        kselfSet: list, default=None
            Manual KPOINTS settings. If kselfSet is a list, the length of the kselfSet should be 3.
        kmeshrv: float, default=0
            The k-mesh accuracy level.
        gamma: bool, default=False
            Whether to use the gamma-centered k-point.
        postype: str, default=None
            For generating the KPOINTS file for band structure calculation.
            If postype is not None and inputOpt includes 'band', the line-mode k-path is generated.
            postype should be 'H' (hexagonal) or 'R' (rectangluar).
        angle: float, default=120
            The angle of the hexagonal lattice.
        is2D: bool, default=True
            Whether the structure is a 2D material.
            If is2D is True, the K-mesh in the z-direction is set to 1.
        potpath: str, default=None
            The path of the POTCAR file.
        is_print: bool, default=False
            Whether to print the path where the input files are written.
        taskname: str, default=None
            The name of the task, which will be used as the suffix of the input file name.
        """
        if stfile is None:
            stfile = []
            for stf_ in os.listdir(posdir):
                stfile.append(os.path.join(posdir, stf_))
        self.StFile = stfile
        self.posdir = posdir
        self.workdir = workdir
        self.overwrite = overwrite
        self.cal_mode = inputOpt
        self.user_incar_settings = user_incar_settings
        self.ini_magmom = ini_magmom
        self.magmom_setting = magmom_setting
        self.nonmagmom = nonmagmom
        self.kselfSet = kselfSet
        self.kmeshrv = kmeshrv
        self.gamma = gamma
        self.postype = postype
        self.angle = angle if angle is not None else 120
        self.is2D = is2D
        self.potpath = potpath
        self.is_print = is_print
        self.taskname = taskname.lower() if taskname is not None else None
        # if self.taskname is not None:
        #     self.multitask = True
        # else:
        #     self.multitask = False

        self.runfile = os.path.join(self.workdir, "run.sh")
        self.StObj = None

    def __write_input(self, workdir_):
        if self.postype is not None and 'band' in self.cal_mode:
            BandKPOINTSet(workdir_, postype=self.postype, angle=self.angle).writefile()
        else:
            MeshKPOINTSet(workdir_, kselfSet=self.kselfSet, StObj=self.StObj, kmeshrv=self.kmeshrv, gamma=self.gamma,
                          is2D=self.is2D).writefile()
        POTCARSet(workdir_, elements=self.StObj.elements, potpath=self.potpath).writefile()
        if self.cal_mode.lower() == 'vepu':
            INCARSetVePu(os.path.join(workdir_, "POSCAR"), os.path.join(workdir_, "POTCAR"), workdir_).writefile()
        elif self.cal_mode.lower() == 'pu':
            INCARSetVePu(os.path.join(workdir_, "POSCAR"), os.path.join(workdir_, "POTCAR"), workdir_, skipVe=True).writefile()
        else:
            INCARSet(workdir_, stfile=self.StObj, inputOpt=self.cal_mode, user_incar_settings=self.user_incar_settings,
                     ini_magmom=self.ini_magmom, magmom_setting=self.magmom_setting, nonmagmom=self.nonmagmom).writefile()
        if self.is_print:
            print(f"The input files are written in the directory: {workdir_}")

        if self.taskname is not None:
            for inpfile in ["INCAR", "KPOINTS"]:
                if os.path.exists(os.path.join(workdir_, inpfile + "-" + self.taskname)):
                    os.remove(os.path.join(workdir_, inpfile + "-" + self.taskname))
                os.rename(os.path.join(workdir_, inpfile),
                          os.path.join(workdir_, inpfile + "-" + self.taskname))

    def write_single_input(self, stf=None):
        self.StObj = StructureFromFile(stf)
        workdir_ = os.path.join(self.workdir, self.StObj.uid)
        if not self.overwrite:
            while os.path.isdir(workdir_):
                self.StObj.uid += '-n'
                workdir_ = os.path.join(self.workdir, self.StObj.uid)
        else:
            if not os.path.isdir(workdir_):
                # shutil.rmtree(workdir_)
                os.makedirs(workdir_)
        copyfile(stf, f'{workdir_}/POSCAR') if stf != f'{workdir_}/POSCAR' else None
        # self.workdir_list.append(workdir_)
        self.__write_input(workdir_)

    def write_input(self):
        if isinstance(self.StFile, list):
            for stf_ in self.StFile:
                self.write_single_input(stf_)
        elif isinstance(self.StFile, str):
            self.write_single_input(self.StFile)
        else:
            raise ValueError("Please provide the valid POSCAR file.")

    def write_input_multi_pos(self):
        """
        Write input files for multiple POSCAR files (POSCAR files saved in multi-level directories).
        Parameters {stfile}, {workdir} and {overwrite} are not used in this function. All operations are based on {posdir}.
        """
        assert os.path.isdir(self.posdir), f"Direcotry {self.posdir} does not exist."
        self.runfile = os.path.join(self.posdir, 'run.sh')

        std = os.listdir(self.posdir)
        work_list = []
        stf_list = []
        for std_ in std:
            base_dir = os.path.join(self.posdir, std_)
            for root, dirs, files in os.walk(base_dir):
                if 'POSCAR' in files:
                    stf_list.append(os.path.join(root, 'POSCAR'))
                    work_list.append(base_dir)
                    break

        for stf_i, stf in enumerate(stf_list):
            self.StObj = StructureFromFile(stf)
            self.__write_input(work_list[stf_i])

    def get_incar_set(self):
        """
        Get the optional INCAR parameters in INCARPara.json.
        """
        opts = list(INCAR_para.keys())
        opts.remove("LISTPara")
        print(f"The available calculation modes are as follows: {opts}")
        return opts

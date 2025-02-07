import itertools

# Common valence states of elements
import os

valance_state = {
    "H": [1],
    "C": [-4],
    "N": [-3],
    "O": [-2],
    "F": [-1],
    "Cl": [-1],
    "Br": [-1],
    "I": [-1],
    "P": [5, -3],
    "S": [-2],

    "Sc": [3],
    "Ti": [3, 4],
    "V": [2, 3, 4, 5],
    "Cr": [3, 2, 6],
    "Mn": [2, 3, 4, 6, 7],
    "Fe": [3, 2],
    "Co": [2, 3],
    "Ni": [2],
    "Cu": [2, 1],
    "Zn": [2],
    "Ga": [3],
    "Ge": [4, 2],

    "Y": [3],
    "Zr": [4],
    "Nb": [5],
    "Mo": [4, 6],
    "Tc": [4, 7],
    "Ru": [3, 4],
    "Rh": [3],
    "Pd": [2, 4],
    "Ag": [1],
    "Cd": [2],
    "In": [3],
    "Sn": [4, 2],
    "Sb": [3, 5],

    "Hf": [4],
    "Ta": [5],
    "W": [4, 6],
    "Re": [4, 7],
    "Os": [4, 6, 8],
    "Ir": [3, 4],
    "Pt": [2, 4],
    "Au": [3, 1],
    "Hg": [1, 2],
    "Tl": [1, 3],
    "Pb": [4, 3, 2],
    "Bi": [3, 5],
    "La": [3],
    "Ce": [3, 4],
    "Pr": [3, 4],
    "Nd": [3],
    "Lu": [3],
}

# The period of the elements for setting U value
period_info = {
    "Sc": "3d",
    "Ti": "3d",
    "V": "3d",
    "Cr": "3d",
    "Mn": "3d",
    "Fe": "3d",
    "Co": "3d",
    "Ni": "3d",
    "Cu": "3d",
    "Zn": "3d",
    "Ga": "3d",
    "Ge": "3d",

    "Y": "4d",
    "Zr": "4d",
    "Nb": "4d",
    "Mo": "4d",
    "Tc": "4d",
    "Ru": "4d",
    "Rh": "4d",
    "Pd": "4d",
    "Ag": "4d",
    "Cd": "4d",
    "In": "4d",
    "Sn": "4d",
    "Sb": "4d",

    "Hf": "5d",
    "Ta": "5d",
    "W": "5d",
    "Re": "5d",
    "Os": "5d",
    "Ir": "5d",
    "Pt": "5d",
    "Au": "5d",
    "Hg": "5d",
    "Tl": "5d",
    "Pb": "5d",
    "Bi": "5d",
    "La": "5d",
    "Ce": "5d",
    "Pr": "5d",
    "Nd": "5d",
    "Lu": "5d",
}

# The U value for the d elements
U_info = {"3d": [2, 4.0, 1.0], "4d": [2, 2.0, 1.0], "5d": [2, 1.0, 1.0], 'other': [-1, 0.0, 0.0]}

# The basic INCAR setting
BAISC_PARA = """PREC   =    Normal
ENCUT  =    500
LREAL  =    Auto
ISTART =    0
ICHARG =    2
NWRITE =    1
NPAR   =    4
ISPIN  =    2
ISYM   =    0
GGA    =    PE
VOSKOWN=    1
IVDW   =     12
IDIPOL =     3
DIPOL  =     0.50000000 0.50000000 0.50000000
LDIPOL =     .TRUE.
NELM   =    200
NELMIN =    6
EDIFF  =    1E-06
IALGO   =    38
ISIF   =    2
EDIFFG =    -0.01
NSW    =    300
IBRION =    2
POTIM  =    0.5
AMIX    = 0.2
BMIX    = 0.0001
AMIX_MAG= 0.8
BMIX_MAG= 0.0001
LMAXMIX =  4
ISMEAR =    0
SIGMA  =   0.05
LORBIT =   11
LWAVE  =    .FALSE.
LCHARG =    .FALSE.
"""


class INCARSetVePu:
    """
    Calculate the valence of the system and set the correct number of electrons in INCAR based on commonly used valences.
    In addition, the U-value is set according to the period in which the element is located.
    Please note that POTCAR and POSCAR are required.

    poscar: str
        The path of the POSCAR file.
    potcar: str
        The path of the POTCAR file.
    save_path: str
        The path to save the INCAR file.
    skipVe: bool
        Whether to skip the calculation of the valence of the system, i.e., only set the U value.
    """
    def __init__(self, poscar, potcar, save_path=None, skipVe=False):
        self.poscar = poscar
        self.potcar = potcar
        if save_path is None:
            self.save_path = os.getcwd()
        else:
            self.save_path = save_path

        self.elem, self.num = self.poscar_parse()
        self.pot_val = self.potcar_parse()
        self.skipVe = skipVe

    def poscar_parse(self):
        assert os.path.exists(self.poscar), "The POSCAR file does not exist!"
        with open(self.poscar, "r") as f:
            inf = f.readlines()
        return inf[5].split(), inf[6].split()

    def potcar_parse(self):
        assert os.path.exists(self.potcar), "The POTCAR file does not exist!"
        pot_val = []
        with open(self.potcar, "r") as f:
            potcar = f.readlines()
            for line in potcar:
                if "ZVAL   =" in line:
                    ptv = line.split()[5]
                    pot_val.append(float(ptv))
        return pot_val

    def get_pot_et(self):
        """
        According to the POTCAR file, get the total electron number of the system.
        """
        pot_et = 0
        for i, e_ in enumerate(self.pot_val):
            pot_et += e_ * int(self.num[i])
        return pot_et

    def get_chem_et(self):
        """
        According to the POSCAR file, get the total electron number of the system.
        """
        chem_e = []
        for elem_ in self.elem:
            if elem_ in valance_state:
                chem_e.append(valance_state[elem_])

        e_e = [list(x) for x in itertools.product(*chem_e)]
        chem_et = []
        for e_ in e_e:
            et_ = 0
            for i, e__ in enumerate(e_):
                et_ += e__ * int(self.num[i])
            chem_et.append(et_)
        return chem_et

    def get_LDAU_info(self):
        LDAU_info = {"LDAUL": [], "LDAUU": [], "LDAUJ": []}
        for elem_ in self.elem:
            if elem_ in period_info:
                LDAU_info["LDAUL"].append(U_info[period_info[elem_]][0])
                LDAU_info["LDAUU"].append(U_info[period_info[elem_]][1])
                LDAU_info["LDAUJ"].append(U_info[period_info[elem_]][2])
            else:
                LDAU_info["LDAUL"].append(U_info["other"][0])
                LDAU_info["LDAUU"].append(U_info["other"][1])
                LDAU_info["LDAUJ"].append(U_info["other"][2])
        return LDAU_info

    def writefile(self):
        pot_et = self.get_pot_et()
        if not self.skipVe:
            chem_et = self.get_chem_et()
            nele = int(pot_et - chem_et[0])
        LDAU_info = self.get_LDAU_info()

        savefile = os.path.join(self.save_path, "INCAR")
        with open(savefile, 'w', newline='\n') as wf:
            wf.write(BAISC_PARA)
            if not self.skipVe:
                wf.write("NELECT = " + str(nele) + "\n")
            wf.write("LDAU = .TRUE.\n")
            wf.write("LDAUTYPE = 2\n")
            wf.write("LDAUL = ")
            for i in range(len(self.elem)):
                wf.write(str(LDAU_info["LDAUL"][i]) + " ")
            wf.write("\n")
            wf.write("LDAUU = ")
            for i in range(len(self.elem)):
                wf.write(str(LDAU_info["LDAUU"][i]) + " ")
            wf.write("\n")
            wf.write("LDAUJ = ")
            for i in range(len(self.elem)):
                wf.write(str(LDAU_info["LDAUJ"][i]) + " ")
            wf.write("\n")
        if not self.skipVe:
            for chem_et_ in chem_et[1:]:
                nele = int(pot_et - chem_et_)
                savefile = os.path.join(self.save_path, "INCAR_" + str(chem_et_))
                with open(savefile, 'w', newline='\n') as wf:
                    wf.write(BAISC_PARA)
                    wf.write("NELECT = " + str(nele) + "\n")
                    wf.write("LDAU = .TRUE.\n")
                    wf.write("LDAUTYPE = 2\n")
                    wf.write("LDAUL = ")
                    for i in range(len(self.elem)):
                        wf.write(str(LDAU_info["LDAUL"][i]) + " ")
                    wf.write("\n")
                    wf.write("LDAUU = ")
                    for i in range(len(self.elem)):
                        wf.write(str(LDAU_info["LDAUU"][i]) + " ")
                    wf.write("\n")
                    wf.write("LDAUJ = ")
                    for i in range(len(self.elem)):
                        wf.write(str(LDAU_info["LDAUJ"][i]) + " ")
                    wf.write("\n")



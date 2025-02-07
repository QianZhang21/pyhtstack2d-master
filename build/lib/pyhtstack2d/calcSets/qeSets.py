import os
import json

import numpy as np

Elem_json = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'periodic_table.json')
with open(Elem_json, "r", encoding='utf-8') as f:
    Element_dic = json.load(f)
Bohr = 0.52917721067


def para_add_spin(para):
    para_spin = para.copy()
    para_spin["system"]["nspin"] = 2
    para_spin["system"]["starting_magnetization(1)"] = 0.5
    return para_spin


QEBlock_json = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QEBlock.json')
with open(QEBlock_json, "r", encoding='utf-8') as f:
    QEBlock = json.load(f)

QeBasicParamRelax_json = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QEParaRelax.json')
with open(QeBasicParamRelax_json, "r", encoding='utf-8') as f:
    QeBasicParamRelax = json.load(f)
QeBasicParamRelaxSpin = para_add_spin(QeBasicParamRelax)

QeBasicParamScf_json = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QEParaScf.json')
with open(QeBasicParamScf_json, "r", encoding='utf-8') as f:
    QeBasicParamScf = json.load(f)
QeBasicParamScfSpin = para_add_spin(QeBasicParamScf)

QeBasicParamMD_json = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'QEParaMD.json')
with open(QeBasicParamMD_json, "r", encoding='utf-8') as f:
    QeBasicParamMD = json.load(f)
QeBasicParamMDSpin = para_add_spin(QeBasicParamMD)


class InputWriter4qe:
    def __init__(self, stfile=None, posdir="POSCAR_dir", workdir="work", inputOpt='scf', spinOpen=False,
                 user_settings=None, block_setting=None, overwrite=True, taskname=None):
        """
        Generate input files for Quantum Espresso (QE) according to the POSCAR file.

        stfile: str or list
            The POSCAR file name or the list of POSCAR file names.
            When stfile is None, it will search the POSCAR files in the posdir.
        posdir: str
            The directory where the POSCAR files are stored.
        workdir: str
            The directory where the input files are stored.
        inputOpt: str
            The type of input file to be generated. Options include 'scf', 'relax' and 'md'.
        spinOpen: bool
            Whether the input file is for spin-polarized calculation. Default is False.
        user_settings: pos_dict
            The user-defined settings for the input file.
        block_setting: pos_dict
            The settings for writing the block tiles and ending tags.
            The basic content of block_setting is as follows:
            {"control": {"blockname": "&CONTROL", "end": "/"}, ...}
        overwrite: bool
            Whether to overwrite the existing files. Default is True.
        taskname: str
            The name of the task, which will be used as the prefix of the input file name.
        """
        if stfile is None:
            stfile = []
            for stf_ in os.listdir(posdir):
                stfile.append(os.path.join(posdir, stf_))
        self.StFile = stfile
        self.posdir = posdir
        self.workdir = workdir
        self.inputOpt = inputOpt.lower()
        self.SpinOpen = spinOpen
        self.overwrite = overwrite
        self.ini_para()
        self.QEBlock = QEBlock.copy()

        if user_settings:
            for key, value in user_settings.items():
                key_ = key.lower()
                if key_ in self.BasicParam:
                    self.BasicParam[key_].update(value)
                else:
                    self.BasicParam[key_] = value

        if block_setting:
            for key, value in block_setting.items():
                key_ = key.lower()
                if key_ in self.QEBlock:
                    self.QEBlock[key_].update(value)
                else:
                    self.QEBlock[key_] = value

        self.runfile = os.path.join(self.workdir, "run.sh")
        self.taskname = taskname.lower() if taskname is not None else None

    def ini_para(self):
        if self.SpinOpen:
            if self.inputOpt == 'scf':
                self.BasicParam = QeBasicParamScfSpin
            elif self.inputOpt == 'relax':
                self.BasicParam = QeBasicParamRelaxSpin
            elif self.inputOpt == 'md':
                self.BasicParam = QeBasicParamMDSpin
                pass
            else:
                raise ValueError("The inputOpt should be 'scf', 'relax' or 'md'. \n For example, the basic parameters "
                                 f"for 'scf' calculation are {QeBasicParamScfSpin}" 
                                 "\n Other options are not supported. Please use user_settings to define specific "
                                 "parameters based on basic parameters.")
        else:
            if self.inputOpt == 'scf':
                self.BasicParam = QeBasicParamScf
            elif self.inputOpt == 'relax':
                self.BasicParam = QeBasicParamRelax
            elif self.inputOpt == 'md':
                self.BasicParam = QeBasicParamMD
            else:
                raise ValueError("The inputOpt should be 'scf', 'relax' or 'md'. \n For example, the basic parameters "
                                 f"for 'scf' calculation are {QeBasicParamScf}" 
                                 "\n Other options are not supported. Please use user_settings to define specific "
                                 "parameters based on basic parameters.")

    def get_para(self):
        print(self.BasicParam)

    def write_single_input(self, stf):
        assert os.path.isfile(stf), f"File {stf} does not exist."
        uid = os.path.basename(stf).split('-')[0]
        workdir_ = os.path.join(self.workdir, uid)
        if not self.overwrite:
            while os.path.isdir(workdir_):
                uid += '-n'
                workdir_ = os.path.join(self.workdir, uid)
        else:
            if not os.path.isdir(workdir_):
                os.makedirs(workdir_)

        with open(stf, "r") as f:
            lines = f.readlines()
        Atom_Type = lines[5].split()
        Atom_Types_Num = len(Atom_Type)

        Atom_Num = lines[6].split()
        Atom_Num_total = 0
        for num in Atom_Num:
            Atom_Num_total = Atom_Num_total + eval(num)

        self.BasicParam["system"]["nat"] = Atom_Num_total
        self.BasicParam["system"]["ntyp"] = Atom_Types_Num
        Atom_Type_list = [item for count, item in zip(Atom_Num, Atom_Type) for _ in range(int(count))]

        cell_factor = float(lines[1].split()[0])
        a = list(map(float, lines[2].split()))
        b = list(map(float, lines[3].split()))
        c = list(map(float, lines[4].split()))

        Super_Cell = cell_factor * np.array([a, b, c])
        Atom_Positions = np.array(
            [list(map(float, line.split()[0:3])) for line in lines[8:8 + Atom_Num_total]])
        if "atomic_positions" not in self.BasicParam:
            self.BasicParam["atomic_positions"] = {}
        pos_type = lines[7].split()[0].lower()
        if pos_type.startswith('d'):
            self.BasicParam["atomic_positions"]["unit"] = "{crystal}"
        elif pos_type.startswith('c'):
            self.BasicParam["atomic_positions"]["unit"] = "{angstrom}"
        else:
            raise ValueError("The atomic position unit is not recognized.")

        with open(os.path.join(workdir_, 'QE.in'), "w") as f:
            for key, value in self.BasicParam.items():
                if key == "atomic_positions":
                    f.write(f"{self.QEBlock[key]['blockname']} {self.BasicParam['atomic_positions']['unit']}\n")
                    for i in range(Atom_Num_total):
                        f.write(f" {Atom_Type_list[i]:<5} {Atom_Positions[i][0]:>15.10f} "
                                f"{Atom_Positions[i][1]:>15.10f} {Atom_Positions[i][2]:>15.10f}\n")
                    f.write(f"{self.QEBlock[key]['end']}")
                elif key == "cell_parameters":
                    f.write(f"{self.QEBlock[key]['blockname']} {self.BasicParam['cell_parameters']['unit']}\n")
                    for i in range(3):
                        f.write(f" {Super_Cell[i][0]:>15.10f} {Super_Cell[i][1]:>15.10f} "
                                f"{Super_Cell[i][2]:>15.10f}\n")
                    f.write(f"{self.QEBlock[key]['end']}")
                elif key == "kpoints":
                    f.write(f"{self.QEBlock[key]['blockname']} {self.BasicParam['kpoints']['ktype']}\n")
                    if 'ktype' in self.BasicParam[key] and self.BasicParam[key]["ktype"] == "{gamma}":
                        pass
                    else:
                        f.write(f"{self.BasicParam[key]['kpoints']}\n")
                    f.write(f"{self.QEBlock[key]['end']}")
                elif key == "pseudo_potential":
                    f.write(f"{self.QEBlock[key]['blockname']}\n")
                    for element in Atom_Type:
                        f.write(f" {element} {Element_dic[element]['atomic_mass']} "
                                f"{element}{self.BasicParam['pseudo_potential']['type']}\n")
                    f.write(f"{self.QEBlock[key]['end']}")
                else:
                    f.write(f"{self.QEBlock[key]['blockname']}\n")
                    for key_, value_ in value.items():
                        f.write(f" {key_} = {value_}\n")
                    f.write(f"{self.QEBlock[key]['end']}\n")
                f.write("\n")

        if self.taskname is not None:
            if os.path.exists(os.path.join(workdir_, self.taskname + "-" + 'QE.in')):
                os.remove(os.path.join(workdir_, self.taskname + "-" + 'QE.in'))
            os.rename(os.path.join(workdir_, 'QE.in'), os.path.join(workdir_, self.taskname + "-" + 'QE.in'))

    def write_input(self):
        if isinstance(self.StFile, list):
            for stf_ in self.StFile:
                self.write_single_input(stf_)
        elif isinstance(self.StFile, str):
            self.write_single_input(self.StFile)
        else:
            raise ValueError("Please provide the valid POSCAR file.")

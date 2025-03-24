import os
import shutil
import json

import numpy as np


Elem_json = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'periodic_table.json')
with open(Elem_json, "r", encoding='utf-8') as f:
    Element_dic = json.load(f)

POT_PATH = os.path.join(os.path.expanduser("~"), ".config", ".PyHTStack2D.json")

yaml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PWDFTconfig.yaml')
Bohr = 0.529177210905667


class InputWriter4pwdft:
    def __init__(self, stfile=None, posdir="POSCAR_dir", workdir="work", user_settings=None, potpath=None,
                 pseudo_suffix="_ONCV_PBE-1.0.upf", overwrite=True, taskname=None):
        """
        Generate input files for PWDFT according to the POSCAR file.

        stfile: str or list
            The POSCAR file name or the list of POSCAR file names.
            When stfile is None, it will search the POSCAR files in the posdir.
        posdir: str
            The directory where the POSCAR files are stored.
        workdir: str
            The directory where the input files are stored.
        user_settings: str
            The path of the user-defined settings file.
            If user_settings is None, the default settings will be used.
        potpath: str
            The path of the pseudo potential file.
        pseudo_suffix: str
            The suffix of the pseudo potential file name.
        overwrite: bool
            Whether to overwrite the existing files. Default is True.
        taskname: str
            The name of the task, which will be used as the prefix of the input file name.
        """
        assert taskname is None, "Multi-tasking with PWDFT is not currently supported."
        if stfile is None:
            stfile = []
            for stf_ in os.listdir(posdir):
                stfile.append(os.path.join(posdir, stf_))
        self.StFile = stfile
        self.posdir = posdir
        self.workdir = workdir
        self.overwrite = overwrite

        if user_settings is None:
            user_settings = yaml_path
        with open(user_settings, 'r') as configf:
            self.user_settings = configf.readlines()

        if potpath is None:
            assert os.path.exists(POT_PATH), "Please check if '~/.config/.PyHTStack2D.json' exists " \
                                             "and make sure the POTCAR path is added."
            with open(POT_PATH, 'r') as f:
                potpath = json.load(f)['POTCAR_PATH']
        self.potpath = potpath
        self.pseudo_suffix = pseudo_suffix
        self.runfile = os.path.join(self.workdir, "run.sh")


    def write_single_input(self, stf):
        assert os.path.isfile(stf), f"File {stf} does not exist."
        uid = ''.join(os.path.basename(stf).split('-')[:-1])
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

        Atom_Type_element_num = []
        for Atom_Type_i in Atom_Type:
            Atom_Type_element_num.append(str(Element_dic[Atom_Type_i]["atomic_no"]))

        Atom_Num = lines[6].split()
        Atom_Num_total = 0
        for num in Atom_Num:
            Atom_Num_total = Atom_Num_total + eval(num)

        cell_factor = float(lines[1].split()[0])
        a = list(map(float, lines[2].split()))
        b = list(map(float, lines[3].split()))
        c = list(map(float, lines[4].split()))

        Super_Cell = cell_factor * np.array([a, b, c])
        Atom_Positions = np.array(
            [list(map(float, line.split()[0:3])) for line in lines[8:8 + Atom_Num_total]])
        pos_type = lines[7].split()[0].lower()
        if pos_type.startswith("c"):
            Atom_Positions = Atom_Positions @ np.linalg.inv(Super_Cell)

        for atom in Atom_Type:
            if not os.path.isfile(os.path.join(self.potpath, atom + self.pseudo_suffix)):
                raise FileNotFoundError(f"Pseudo potential file {atom + self.pseudo_suffix} does not exist.")
            if os.path.exists(os.path.join(workdir_, atom + self.pseudo_suffix)):
                pass
            else:
                shutil.copyfile(os.path.join(self.potpath, atom + self.pseudo_suffix),
                                os.path.join(workdir_, atom + self.pseudo_suffix))

        with open(os.path.join(workdir_, 'config.yaml'), "w") as f:
            f.write('UPF_File:\n')
            for atom in Atom_Type:
                f.write(f'  -  {atom}{self.pseudo_suffix}\n')
            f.write('\n')
            for line in self.user_settings:
                f.write(line)
            f.write('\n')
            f.write(f'Atom_Types_Num:        {Atom_Types_Num}\n')
            f.write('Atom_Type:             [ ' + ', '.join(Atom_Type_element_num) + ' ]\n')
            f.write('Super_Cell:\n')
            f.write(f'  -  [ {Super_Cell[0][0]:.10f}, {Super_Cell[0][1]:.10f}, {Super_Cell[0][2]:.10f} ]\n')
            f.write(f'  -  [ {Super_Cell[1][0]:.10f}, {Super_Cell[1][1]:.10f}, {Super_Cell[1][2]:.10f} ]\n')
            f.write(f'  -  [ {Super_Cell[2][0]:.10f}, {Super_Cell[2][1]:.10f}, {Super_Cell[2][2]:.10f} ]\n')

            f.write('Atom_Num:   [ ' + ', '.join(Atom_Num) + ' ]' + '\n')

            f.write('Atom_Red:\n')
            for i in range(Atom_Num_total):
                f.write(f'  -  [ {Atom_Positions[i][0]:.10f}, {Atom_Positions[i][1]:.10f}, {Atom_Positions[i][2]:.10f} ]\n')

    def write_input(self):
        if isinstance(self.StFile, list):
            for stf_ in self.StFile:
                self.write_single_input(stf_)
        elif isinstance(self.StFile, str):
            self.write_single_input(self.StFile)
        else:
            raise ValueError("Please provide the valid POSCAR file.")


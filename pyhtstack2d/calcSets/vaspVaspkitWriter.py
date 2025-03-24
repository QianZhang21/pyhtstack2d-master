"""
This module is used to write VASP input files for a given POSCAR file using the VASPKIT package.
About the VASPKIT package, please refer to the official website: https://vaspkit.com/
"""
import os
# import shutil
import warnings
from shutil import copyfile


class VaspkitInputWriter:
    def __init__(self, stfile=None, posdir="POSCAR_dir", workdir="work", overwrite=True, inputOpt='ST', incexis=False,
                 kmseshScheme=2, kmeshrv=0.04, bandkpath=None, Hybridbandkpath=None, kpexis=False,
                 genpotcar=103, potexis=False, is_print=False, taskname=None, shoverwrite=False):
        """
        Use the VASPKIT package to write VASP input files for a given POSCAR file. The VASPKIT package is a powerful tool
        for pre- and post-processing VASP calculations.

        Details of the VASPKIT package can be found at: https://vaspkit.com/

        stfile: str or list
            The POSCAR file name or a list of POSCAR file names.
        posdir: str
            The directory where the POSCAR files are located.
            For example, posdir = 'POSCAR_dir' means the POSCAR files are located in the 'POSCAR_dir' directory.
        workdir: str
            The directory where the input files are written.
        overwrite: bool
            Whether to overwrite the existing input files.
        inputOpt: str
            INCAR Options. Consistent with the menu of the VASPKIT tool generation template INCAR Options.
            The VASPKIT common "101" shows the available options, including: ST, SR, MG, SO, D3,
            H6, PU, MD, GW, BS, DC, EL, BD, OP, EC, PC, FD, DT, NE, DM, FQ and LR.
            The specific meaning of each option can be found in the VASPKIT package documentation.
            For example, inputOpt = 'STH6D3' means HSE06-D3 Static-Calcualtion'.
        incexis: bool
            Whether the INCAR file exists.
            If True, skip the INCAR file generation.
        kmseshScheme: int
            K-Mesh Scheme.
            1: Monkhorst-Pack Scheme; 2: Gamma-Centered Scheme; 3: Irreducible K-Points with Gamma Scheme.
        kmeshrv: float
            Kmesh-Resolved Value. Accuracy level: Gamma-Only: 0; Low: 0.06~0.04; Medium: 0.04~0.03; Fine: 0.02-0.01.
        bandkpath: int
            K-Path Options for Band Structure Calculation.
            301: 1D Structure; 302: 2D Structure; 303: 3D Structure.
            The other options are provided in the VASPKIT package documentation.
        Hybridbandkpath: list
            K-Path Options for Hybrid Band Structure Calculation.
            For example, [251, 2, 0.03, 0.04] means generating KPOINTS for the hybrid band structure calculation
            with Gamma-Centered Scheme, and the k-mesh resolved value for SCF calculation is 0.03, the k-mesh resolved
            value along K-Path for Band Calculation is 0.04.
        kpexis: bool
            Whether the KPOINTS file exists.
            If True, skip the KPOINTS file generation.
        genpotcar: int
            Pseudopotential Generation Scheme.
            103: Generate POTCAR File with Default Setting.
            104: Generate POTCAR File with User Specified Potential
        potexis: bool
            Whether the POTCAR file exists.
            If True, skip the POTCAR file generation.
        is_print: bool
            Whether to print the path where the input files are written.
        taskname: str
            The name of the task, which will be used as the suffix of the input file name.
        shoverwrite: bool
            Whether to overwrite the existing mkdir.sh file.
        """
        assert kmseshScheme in [1, 2, 3], "kmseshScheme must be an integer and in [1, 2, 3]."
        assert isinstance(Hybridbandkpath, list) and len(Hybridbandkpath) == 4 if Hybridbandkpath is not None else True, \
            "Hybridbandkpath must be a list with 4 elements. For example, [251, 2, 0.03, 0.04], where 251 is the K-Path" \
            "Options for Hybrid Band Structure Calculation, 2 is the K-Mesh Scheme, 0.03 is the k-mesh resolved value " \
            "for SCF calculation, and 0.04 is the k-mesh resolved value along K-Path for Band Calculation."
        assert isinstance(genpotcar, int) and genpotcar in [103, 104], "genpotcar must be an integer and in [103, 104]."

        if Hybridbandkpath is not None and bandkpath is None:
            bandkpath = 302
            warnings.warn("If Hybridbandkpath is not None and bandkpath is missing, bandkpath is set to 302 (2D structure) for KPATH.in generation.")

        if stfile is None:
            stfile = []
            for stf_ in os.listdir(posdir):
                stfile.append(os.path.join(posdir, stf_))
        self.StFile = stfile
        self.posdir = posdir
        self.workdir = workdir
        self.overwrite = overwrite
        self.inputOpt = inputOpt
        self.incexis = incexis
        self.kmseshScheme = kmseshScheme
        self.krv = kmeshrv
        self.bandkpath = bandkpath
        self.Hybridbandkpath = Hybridbandkpath
        self.kpexis = kpexis
        self.genpotcar = genpotcar
        self.potexis = potexis
        self.is_print = is_print
        self.taskname = taskname.lower() if taskname is not None else None

        # self.workdir_list = []
        self.bashfile = os.path.join(self.workdir, 'mkdir.sh')
        self.runfile = os.path.join(self.workdir, 'run.sh')
        self.shoverwrite = shoverwrite
        if not self.shoverwrite:
            if not os.path.isfile(self.bashfile):
                warnings.warn(f"File {self.bashfile} does not exist. The mkdir.sh file will be generated.")
                self.shoverwrite = True

    def __write_input(self):
        # if os.path.isfile(self.runfile):
        #     os.remove(self.runfile)
        if not os.path.isdir(os.path.dirname(self.bashfile)):
            os.makedirs(os.path.dirname(self.bashfile))

        if self.shoverwrite or not os.path.isfile(self.bashfile):
            with open(self.bashfile, 'w', newline='\n') as f:
                f.write('#!/bin/bash\n\n')
        else:
            with open(self.bashfile, 'a', newline='\n') as f:
                f.write('\n')

        with open(self.bashfile, 'a', newline='\n') as f:
            f.write(f"for file in `ls $1`\n")
            f.write(f"do\n")
            f.write(f"    if [ -d $file ]; then\n")
            f.write(f"        cd $file\n")
            if not self.incexis:
                f.write(f'        echo -e "101\\n{self.inputOpt}" | vaspkit >/dev/null 2>&1\n')

            if not self.potexis:
                f.write(f'        echo -e "{self.genpotcar}" | vaspkit >/dev/null 2>&1\n')

            if not self.kpexis:
                if self.Hybridbandkpath is not None:
                    f.write(f'        echo -e "{self.bandkpath}" | vaspkit >/dev/null 2>&1\n')
                    f.write(
                        f'        echo -e "{self.Hybridbandkpath[0]}\\n{self.Hybridbandkpath[1]}\\n{self.Hybridbandkpath[2]}\\n'
                        f'{self.Hybridbandkpath[3]}" | vaspkit >/dev/null 2>&1\n')
                elif self.bandkpath is not None:
                    f.write(f'        echo -e "{self.bandkpath}" | vaspkit >/dev/null 2>&1\n')
                    f.write('        cp KPATH.in KPOINTS\n')
                else:
                    f.write(f'        echo -e "102\\n{self.kmseshScheme}\\n{self.krv}" | vaspkit >/dev/null 2>&1\n')

            if self.taskname is not None:
                if not self.incexis:
                    f.write(f'        mv INCAR INCAR-{self.taskname}\n')
                if not self.kpexis:
                    f.write(f'        mv KPOINTS KPOINTS-{self.taskname}\n')
            f.write('        cd ../\n')
            f.write('    fi\n')
            f.write('done\n')

    def write_single_input(self, stf=None):
        assert os.path.isfile(stf), f"File {stf} does not exist."
        uid = ''.join(os.path.basename(stf).split('-')[:-1])
        workdir_ = os.path.join(self.workdir, uid)
        if not self.overwrite:
            while os.path.isdir(workdir_):
                uid += '-n'
                workdir_ = os.path.join(self.workdir, uid)
        else:
            if not os.path.isdir(workdir_):
                # shutil.rmtree(workdir_)
                os.makedirs(workdir_)
        copyfile(stf, f'{workdir_}/POSCAR') if stf != f'{workdir_}/POSCAR' else None
        # self.workdir_list.append(workdir_)
        if self.is_print:
            print(f'Input files are written in {workdir_}')

    def write_input(self):
        if isinstance(self.StFile, list):
            for stf_ in self.StFile:
                self.write_single_input(stf_)
        elif isinstance(self.StFile, str):
            self.write_single_input(self.StFile)
        else:
            raise ValueError("Please provide the valid POSCAR file.")

        self.__write_input()

    def write_input_multi_pos(self):
        """
        Write input files for multiple POSCAR files (POSCAR files saved in multi-level directories).
        Parameters {stfile}, {workdir} and {overwrite} are not used in this function. All operations are based on {posdir}.
        """
        assert os.path.isdir(self.posdir), f"Direcotry {self.posdir} does not exist."
        self.bashfile = os.path.join(self.posdir, 'mkdir.sh')
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
            copyfile(stf, f'{work_list[stf_i]}/POSCAR') if stf != f'{work_list[stf_i]}/POSCAR' else None

        self.__write_input()




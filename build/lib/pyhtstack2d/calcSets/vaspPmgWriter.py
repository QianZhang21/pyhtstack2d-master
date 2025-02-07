"""
This module, which depends on the pymatgen package, provides a class for writing VASP input files for a given POSCAR file.
About the pymatgen package, please refer to the official website: https://pymatgen.org/
"""
import os

from pymatgen.io.vasp.sets import *

setModes = {"MITMD": MITMDSet, "MITNEB": MITNEBSet, "MITRelax": MITRelaxSet, "MPAbsor": MPAbsorptionSet,
            "MPHSEBS": MPHSEBSSet, "MPHSERelax": MPHSERelaxSet, "MPMD": MPMDSet, "MPMetalRelax": MPMetalRelaxSet,
            "MPNMR": MPNMRSet, "MPNonSCF": MPNonSCFSet, "MPRelax": MPRelaxSet, "MPSOC": MPSOCSet,
            "MPScanRelax": MPScanRelaxSet, "MPScanStatic": MPScanStaticSet, "MPStatic": MPStaticSet,
            "MSONable": MSONable, "MVLElastic": MVLElasticSet, "MVLGB": MVLGBSet, "MVLGW": MVLGWSet,
            "MVLNPTMD": MVLNPTMDSet, "MVLRelax52": MVLRelax52Set, "MVLScanRelax": MVLScanRelaxSet,
            "MVLSlab": MVLSlabSet, "MatPESStatic": MatPESStaticSet}

pmappath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pmap.json')
pmap = loadfn(pmappath)


class PmgInputWriter:
    def __init__(self, stfile=None, posdir="POSCAR_dir", workdir="work", inputOpt="MPRelax", overwrite=True, user_incar_settings=None,
                 is_print=False, taskname=None):
        """
        Write VASP input files for a given POSCAR file. The pymatgen package provides a set of classes to generate VASP
        input files. Note that the location of the pseudopotential file has been pre-defined, i.e. the pmgrc.yaml file
        contains the path where the pseudopotential file is located.

        About the pymatgen package, please refer to the official website: https://pymatgen.org/

        stfile: str
            Path to the POSCAR file.
        posdir: str
            The directory where the POSCAR files are located.
            For example, posdir = 'POSCAR_dir' means the POSCAR files are located in the 'POSCAR_dir' directory.
        workdir: str
            The directory where the VASP input files are written.
        inputOpt: str
            The style of the settings in the VASP input file.
            The available styles are: "MITMD", "MITNEB", "MITRelax", "MPAbsor", "MPHSEBS", "MPHSERelax", "MPMD",
            "MPMetalRelax", "MPNMR", "MPNonSCF", "MPRelax", "MPSOC", "MPScanRelax", "MPScanStatic", "MPStatic",
            "MSONable", "MVLElastic", "MVLGB", "MVLGW", "MVLNPTMD", "MVLRelax52", "MVLScanRelax", "MVLSlab",
            "MatPESStatic".
        overwrite: bool
            Whether to overwrite the existing input files.
        user_incar_settings: dict
            User-defined INCAR settings.
            Generally, 99% of tasks require specialization of the INCAR file by user_incar_settings.
            For example, user_incar_settings = {"NPAR": 4, "EDIFF": 1e-5, , "IVDW": 12,
            "LDIPOL": "True", "IDIPOL": 3,  "DIPOL": [0.5, 0.5, 0.5]}
        is_print: bool
            Whether to print the path where the input files are written.
        tasklist: str
            The name of the task, which will be used as the suffix of the input file name.
        """
        assert inputOpt in setModes, f"Set mode {inputOpt} is not supported. \n" \
                                     f" The available modes are  {setModes.keys()}"
        if stfile is None:
            stfile = []
            for stf_ in os.listdir(posdir):
                stfile.append(os.path.join(posdir, stf_))
        self.StFile = stfile
        self.posdir = posdir
        self.workdir = workdir
        self.overwrite = overwrite
        self.user_incar_settings = user_incar_settings
        self.is_print = is_print
        self.taskname = taskname.lower() if taskname is not None else None
        if self.taskname is not None:
            self.multitask = True
        else:
            self.multitask = False
        self.pmgSetMode = setModes[inputOpt]
        self.SetModesDict = list(setModes.keys())

        # self.workdir_list = []
        self.runfile = os.path.join(self.workdir, 'run.sh')

    def __write_input(self, struct, workdir_):
        input_file = self.pmgSetMode(struct, user_incar_settings=self.user_incar_settings)
        inc = input_file.incar
        kpo = input_file.kpoints
        pos = Poscar(structure=input_file.structure)
        symbols = [pmap[el] for el in pos.site_symbols]
        try:
            pot = Potcar(symbols=symbols)
            viw = VaspInput(poscar=pos, incar=inc, potcar=pot, kpoints=kpo)
            viw.write_input(workdir_)

            if self.multitask:
                for inpfile in ["INCAR", "KPOINTS"]:
                    if os.path.exists(os.path.join(workdir_, inpfile + "-" + self.taskname)):
                        os.remove(os.path.join(workdir_, inpfile + "-" + self.taskname))
                    os.rename(os.path.join(workdir_, inpfile),
                              os.path.join(workdir_, inpfile + "-" + self.taskname))
        except:
            print("Error: The pseudopotential file is not found. Please check the pmgrc.yaml file. \n"
                  "You may use the commands 'pmg config -p /path/to/existing/pseudopotential "
                  "/path/to/your/pseudopotential' \n"
                  "and 'pmg config --add PMG_VASP_PSP_DIR /path/to/your/pseudopotential' "
                  "to add the path to the pseudopotential file.")

    def write_single_input(self, stf=None):
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
        # self.workdir_list.append(workdir_)
        try:
            struct = Structure.from_file(stf)
            try:
                struct.to(os.path.join(workdir_, uid + '.vasp'), 'poscar')
            except:
                struct.to('poscar', os.path.join(workdir_, uid + '.vasp'))

            self.__write_input(struct, workdir_)

            if self.is_print:
                print(f"Input files for {uid} have been written to {workdir_}.")

        except:
            print(f"Error: {stf} is not a valid POSCAR file.")

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
                    work_list.append(base_dir)
                    stf_list.append(os.path.join(root, 'POSCAR'))
                    break

        for stf_i, stf in enumerate(stf_list):
            struct = Structure.from_file(stf)
            self.__write_input(struct, work_list[stf_i])
            os.remove(os.path.join(work_list[stf_i], 'POSCAR'))
            if self.is_print:
                print(f"Input files have been written to {work_list[stf_i]}.")

    def mode_show(self):
        print("The VASP setting modes available in the pymatgen library are:")
        for i_ in range(int(len(self.SetModesDict) / 6)):
            print("  ".join(key for key in self.SetModesDict[i_ * 6:(i_ + 1) * 6]))
        # print("  ".join(key for key in self.SetModesDict[(i_+1)*6:]))
        print(f"Please choose one of the above modes for the VASP input file settings.")

import os
import shutil
import warnings


class ShareInp:
    def __init__(self, workdir, incpath=None, kppath=None, potpath=None, tasklist=None):
        """
        Self-provisioning input file paths for sharing in different calculation directories.

        workdir: str
            The path of the working directory.
        incpath: list
            The path list of the INCAR file, for example ["path/to/INCAR_opt", "path/to/INCAR_scf", "path/to/INCAR_band"].
            The order of the example list is as follows: INCAR for Structural Optimisation, INCAR for Static Calculations, INCAR for Band Calculations.
        kpopath: list
            The path list of the KPOINTS file, for example ["path/to/KPOINTS", "path/to/KPOINTS", "path/to/KPATH.in"].
        potpath: str
            The path of the POTCAR file.
        tasklist: list
            The task name, such as ['relax', 'static', 'band'].
        """
        self.incshare = False
        self.kpshare = False
        self.potshare = False
        self.ntask = len(tasklist) if tasklist is not None else 0

        if incpath is not None:
            self.ntask = len(incpath) if len(incpath) < len(tasklist) else len(tasklist)
            if len(incpath) != len(tasklist):
                warnings.warn(f"The length of the INCAR path list should be equal to the tasklist.\n"
                              f"Only the {tasklist[:self.ntask]} tasks will be created.")
            assert all([os.path.exists(inc) for inc in incpath[:self.ntask]]), "The INCAR file does not exist."
            self.incshare = True
        if kppath is not None:
            self.ntask = len(kppath) if len(kppath) < len(tasklist) else len(tasklist)
            if len(kppath) != len(tasklist):
                warnings.warn(f"The length of the KPOINTS path list should be equal to the tasklist.\n"
                              f"Only the {tasklist[:self.ntask]} tasks will be created.")
            assert all([os.path.exists(kp) for kp in kppath[:self.ntask]]), "The KPOINTS file does not exist."
            self.kpshare = True
        if potpath is not None:
            assert os.path.exists(potpath), "The POTCAR file does not exist."
            self.potshare = True

        self.workdir = workdir
        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        self.INCAR_path = incpath
        self.KPOINTS_path = kppath
        self.POTCAR_path = potpath
        if self.ntask != 0:
            self.tasklist = [task_.lower() for task_ in tasklist]
        else:
            self.tasklist = []

        self.ini_share_file()

    def ini_share_file(self):
        if self.incshare:
            for i, task_ in enumerate(self.tasklist[:self.ntask]):
                if self.INCAR_path[i] != os.path.join(self.workdir, "INCAR-" + task_):
                    shutil.copy(self.INCAR_path[i], os.path.join(self.workdir, "INCAR-" + task_))
        if self.kpshare:
            for i, task_ in enumerate(self.tasklist[:self.ntask]):
                if self.KPOINTS_path[i] != os.path.join(self.workdir, "KPOINTS-" + task_):
                    shutil.copy(self.KPOINTS_path[i], os.path.join(self.workdir, "KPOINTS-" + task_))
        if self.potshare:
            if self.POTCAR_path != os.path.join(self.workdir, "POTCAR"):
                shutil.copy(self.POTCAR_path, os.path.join(self.workdir, "POTCAR"))

    def is_share(self):
        return self.incshare, self.kpshare, self.potshare

import os
import json
import subprocess
import warnings

import numpy as np

from .analysisPROCAR import readPROCAR
from .plotBAND import plotBS


def run_subprocess(command):
    return subprocess.run(command, shell=True, capture_output=True, text=True).stdout.strip()


def get_band_location(bandpath, band_type, index, suffix=' (TO)'):
    command = f"grep 'Location of {band_type}{suffix}' {bandpath}/BAND_GAP | awk '{{print ${index}}}'"
    return float(run_subprocess(command))


def get_band_index(bandpath, index):
    command1 = f"grep 'Highest-Occupied Band' {bandpath}/BAND_GAP | awk '{{print ${index}}}'"
    command2 = f"grep 'Lowest-Unoccupied Band' {bandpath}/BAND_GAP | awk '{{print ${index}}}'"
    return int(run_subprocess(command1)), int(run_subprocess(command2))


def get_vbm_cbm_locations(bandpath, suffix=' (TO)'):
    if suffix:
        rangestar = 5
    else:
        rangestar = 4
    vbmlocation = [get_band_location(bandpath, "VBM", i, suffix) for i in range(rangestar, rangestar + 3)]
    cbmlocation = [get_band_location(bandpath, "CBM", i, suffix) for i in range(rangestar, rangestar + 3)]
    return vbmlocation, cbmlocation


def get_bm_spin(bandpath):
    up_vbm_band_index, up_cbm_band_index = get_band_index(bandpath, 3)
    vbm_band_index, cbm_band_index = get_band_index(bandpath, 5)
    vbmspin = "1" if up_vbm_band_index == vbm_band_index else "-1"
    cbmspin = "1" if up_cbm_band_index == cbm_band_index else "-1"
    return [vbm_band_index - 1], [cbm_band_index - 1], vbmspin, cbmspin


def is_direct_band_gap(vbmlocation, cbmlocation):
    return all(abs(vbm - cbm) < 0.01 for vbm, cbm in zip(vbmlocation, cbmlocation))


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)


class StructureFile:
    def __init__(self, stfile):
        """
        Get the 2D structure  from a POSCAR file.

        stfile: str
            Path to the POSCAR file.
            Default Z-axis is vacuum layer (out-of-plane) direction.
        """
        self.StFile = stfile
        self.lattice = None
        self.a = None
        self.b = None
        self.c = None
        self.abc = None
        self.elements = None
        self.elemnum = None
        self.posType = None
        self.positions = None
        self.indices = None
        self.positions_direct = None
        self.getInfo()
        self.zdiff = np.diff(self.positions[:, 2])

    def getInfo(self):
        """
        Get the information from the structure file, ensuring it exists and parsing it.
        """
        if not os.path.exists(self.StFile):
            raise FileNotFoundError("The POSCAR file does not exist.")

        with open(self.StFile, "r") as f:
            lines = f.readlines()

        self.lattice = np.array([list(map(float, lines[i].split())) for i in range(2, 5)])
        self.elements = list(map(str, lines[5].split()))
        self.elemnum = list(map(int, lines[6].split()))
        self.posType = lines[7].strip().lower()
        positions_ = np.array([list(map(float, line.split()[:3])) for line in lines[8:8 + sum(self.elemnum)]])
        assert self.lattice.shape == (3, 3), "Lattice vectors are not 3x3."
        assert len(self.elements) == len(self.elemnum), "Number of elements and element numbers do not match."

        self.a = np.linalg.norm(self.lattice[0])
        self.b = np.linalg.norm(self.lattice[1])
        self.c = np.linalg.norm(self.lattice[2])
        self.abc = (self.a, self.b, self.c)

        if self.a > self.c or self.b > self.c:
            raise ValueError("Please check the vacuum layer direction: z-axis is not the longest lattice vector.")

        if self.posType.startswith('d'):
            self.positions_direct = positions_.copy()
            self.positions = np.dot(positions_, self.lattice)
        elif self.posType.startswith('c'):
            self.positions = positions_.copy()
            self.positions_direct = np.dot(positions_, np.linalg.inv(self.lattice))
        else:
            raise ValueError("Please check the position type: it is not cartesian or direct.")

        self.indices = np.argsort(self.positions[:, 2])
        self.positions = self.positions[self.indices]
        self.positions_direct = self.positions_direct[self.indices]
        sorted_elements = []
        for elem, num in zip(self.elements, self.elemnum):
            sorted_elements.extend([elem] * num)
        self.elements = [sorted_elements[i] for i in self.indices]


class GetResults:
    def __init__(self, materialfile="material.txt", scf=None, band=None, hybband=None, optical=None, esp=None,
                 directbandgap=False, bandplot=False, banderange=(-3, 3), imag_format="pdf", elempro=True, dim=2,
                 energyunit=1,
                 weightsave=False, nlayer=1, workdir=None, multilevel=None,
                 mag=False, vaspkit=True, pmg=True, infodict=None):
        """
        Extraction of the results after high-throughput calculations by VASP.

        materialfile: str
            Contains paths for information extraction.
        scf: str
            Naming of the file or folder where the opt or SCF calculation results are stored.
        band: str
            Naming of the file or folder where the band structure calculation results are stored.
        hybband: str
            Naming of the file or folder where the hybrid band structure calculation results are stored.
        optical: str
            Naming of the file or folder where the optical properties calculation results are stored.
        esp: str
            Naming of the file or folder where the electrostatic potential calculation results are stored.
        directbandgap: bool
            Whether to get the direct bandgap and the site projected wave function character.
        bandplot: bool
            Whether to plot the band structure, or the layer projected band structure for getInfoBi().
        banderange: tuple
            The energy range for plotting the band structure.
        imag_format: str
            The format for saving the band structure plot.
        elempro: bool
            Whether to extract the element projected band structure / DOS.
        dim: int
            The dimension of the system (1D, 2D or 3D).
        energyunit: int
            The unit of energy for absorption coefficient.
            1: eV; 2: nm; 3: THz; 4: cm-1.
        weightsave: bool
            Whether to save the weight extracted from the PROCAR file.
        nlayer: int
            The number of layers for 2D materials.
        workdir: str
            The working directory where the results are stored.
            If workdir="", also means the current directory.
        multilevel: int
            The number of sub-directories for POSCAR files (stored in a multi-level directory).
            For example, if the POSCAR files are stored in the directory "BiPOSCAR_dir/S24P8Ag4V4_0.00%/AA/cord1", the multilevel is 4.
            If the multilevel is None, it means that each directory (formula) has only one POSCAR file.
        mag: bool
            If True, the magnetic moment is calculated.
        vaspkit: bool
            If True, use VASPKIT to get the information.
        pmg: bool
            If True, use Pymatgen to get the information.
        infodict: json file
            If not None, the information is extracted from the json file and skip the initialization of info_dict by getInfo().
        """
        self.materialfile = materialfile
        self.scf = scf if scf else ""
        self.band = band
        self.hybband = hybband
        self.bandplot = bandplot
        self.directbandgap = directbandgap
        self.banderange = banderange
        self.imag_format = imag_format
        self.elempro = elempro
        self.optical = optical
        self.esp = esp

        self.mag = mag
        assert dim in [1, 2, 3], "The dimension of the system must be 1, 2 or 3."
        self.dim = dim
        assert energyunit in [1, 2, 3, 4], "The unit of energy must be 1, 2, 3 or 4."
        self.energyunit = energyunit
        self.weightsave = True if nlayer == 2 else weightsave
        assert nlayer in [1, 2], "The number of layers must be 1 or 2."
        self.nlayer = nlayer

        self.vaspkit = vaspkit
        self.pmg = pmg
        self.Vasprun, self.BSDOSPlotter, self.BSPlotter = None, None, None
        if self.pmg:
            self.Vasprun = __import__('pymatgen.io.vasp.outputs', fromlist=['Vasprun']).Vasprun
            self.Locpot = __import__('pymatgen.io.vasp.outputs', fromlist=['Locpot']).Locpot
            self.Outcar = __import__('pymatgen.io.vasp.outputs', fromlist=['Outcar']).Outcar
            self.Structure = __import__('pymatgen.core.structure', fromlist=['Structure']).Structure
            self.BSDOSPlotter = __import__('pymatgen.electronic_structure.plotter',
                                           fromlist=['BSDOSPlotter']).BSDOSPlotter
            self.BSPlotter = __import__('pymatgen.electronic_structure.plotter', fromlist=['BSPlotter']).BSPlotter
        elif self.vaspkit:
            print("Please make sure that VASPKIT >= 1.3.1 is installed.")
            if self.directbandgap:
                try:
                    self.Vasprun = __import__('pymatgen.io.vasp.outputs', fromlist=['Vasprun']).Vasprun
                except ImportError:
                    self.Vasprun = None
                    self.directbandgap = False

        if workdir is None:
            self.workdir = os.getcwd()
        elif workdir == "":
            self.workdir = ""
        else:
            assert os.path.exists(workdir), "workdir should be a valid directory."
            self.workdir = workdir

        self.materials = None
        self.infomono_dict = None
        self.infobi_dict = {}
        if multilevel is None or multilevel == 2:
            self.multilevel = None
        else:
            assert isinstance(multilevel, int) and multilevel > 2, "multilevel should be an integer greater than 2."
            self.multilevel = multilevel - 1

        if infodict is not None:
            assert os.path.exists(infodict), "infodict should be a valid file."
            with open(infodict, "r") as f:
                self.info_dict = json.load(f)
        else:
            self.initInfoDict()

    def initInfoDict(self):
        materialfilepath = os.path.join(self.workdir, self.materialfile)
        if not os.path.exists(materialfilepath):
            self.getMaterialfile()
            if not os.path.exists(materialfilepath):
                raise FileNotFoundError(f"{self.materialfile} does not exist.")
        self.materials = self.validmultilevel()

        infosavename = "info"
        infosavenu = 0
        infosavepath = os.path.join(self.workdir, infosavename + ".json")
        while os.path.exists(infosavepath):
            infosavenu += 1
            infosavename = "info" + str(infosavenu)
            infosavepath = os.path.join(self.workdir, infosavename + ".json")
        self.infosavefile = infosavepath
        self.info_dict = {}
        self.getInfo()

    def getMaterialfile(self):
        """
        Get the material file.
        """
        getsh = os.path.join(self.workdir, "getmaterial.sh")
        if self.scf != "" and self.scf.endswith("/"):
            scfdir = self.scf
        elif self.scf != "":
            scfdir = self.scf + "/"
        else:
            scfdir = ""
        with open(getsh, "w", newline="\n") as f:
            f.write("#!/bin/bash\n\n")
            f.write(f"rm -f {self.materialfile}\n\n")
            if self.multilevel is None:
                f.write("for file in $(ls -d */); do\n")
                f.write(f"    if [ -f $file/{scfdir}OUTCAR ] && [ `grep -c \"reached\" $file/{scfdir}OUTCAR` -ne 0 ] && [ `grep -c \"Voluntary\" $file/{scfdir}OUTCAR` -ne 0 ]; then\n")
                f.write(f"        echo $file >> {self.materialfile}\n")
                f.write("    fi\n")
                f.write("done\n")
            else:
                file_path = "$file"
                cd_com = "../" * self.multilevel
                f.write("for file in $(ls -d */); do\n")
                f.write("    cd $file\n")
                for level_i in range(self.multilevel - 1):
                    file_path += f"$file_{level_i}"
                    f.write("    " * level_i + f"    for file_{level_i} in $(ls -d */)\n")
                    f.write("    " * level_i + "    do\n")
                    f.write("    " * level_i + f"        cd $file_{level_i}\n")

                f.write("    " * (self.multilevel - 2) + f"        if [ -f {scfdir}OUTCAR ] && [ `grep -c \"reached\" {scfdir}OUTCAR` -ne 0 ] && [ `grep -c \"Voluntary\" {scfdir}OUTCAR` -ne 0 ]; then\n")
                f.write("    " * (
                            self.multilevel - 2) + f"                echo {file_path} >> {cd_com}{self.materialfile}\n")
                f.write("    " * (self.multilevel - 2) + "        fi\n")
                for level_i in range(self.multilevel - 1):
                    f.write("    " * (self.multilevel - 2 - level_i) + "        cd ..\n")
                    f.write("    " * (self.multilevel - 2 - level_i) + "    done\n")
                f.write("    cd ..\n")
                f.write("done\n")
        os.system("bash " + getsh)

    def validmultilevel(self):
        """
        Read the material file.
        """
        with open(self.materialfile, "r") as f:
            materials = f.readlines()
        material0 = [x for x in materials[0].strip().split("/") if x]

        if len(material0) == 1:
            self.multilevel = None
        else:
            self.multilevel = len(material0)
        return materials

    def getInfo(self):
        def checkpath(path, filename):
            if path is None:
                return None, None
            if not os.path.exists(path):
                print(f"Cannot find the band path {path}.")
                return None, None
            outfile = os.path.join(path, filename)
            if not os.path.exists(outfile):
                print(f"Cannot find the vasprun.xml file in {path}.")
                return None, None
            return path, outfile

        for material in self.materials:
            material = material.strip()
            if material.endswith("/"):
                material = material[:-1]
            posfile = "CONTCAR"
            stfile = os.path.join(self.workdir, material.strip(), posfile)
            if not os.path.exists(stfile):
                posfile = "POSCAR"
                stfile = os.path.join(self.workdir, material.strip(), posfile)
                if not os.path.exists(stfile):
                    raise FileNotFoundError(
                        f"Cannot find the POSCAR/CONTCAR file in {os.path.join(self.workdir, material.strip())}")
            scfoutcarfile = os.path.join(self.workdir, material, self.scf, "OUTCAR")
            scfvasprunxml = os.path.join(self.workdir, material, self.scf, "vasprun.xml")
            scfoszicarfile = os.path.join(self.workdir, material, self.scf, "OSZICAR")
            if not os.path.exists(scfoutcarfile):
                print(f"Cannot find the OUTCAR file in {os.path.join(self.workdir, material, self.scf)}")
                scfoutcarfile = None

            bandpathjoin = os.path.join(self.workdir, material.strip(), self.band) if self.band is not None else None
            bandpath, bandvasprunxml = checkpath(bandpathjoin, "vasprun.xml")

            hybbandpathjoin = os.path.join(self.workdir, material.strip(), self.hybband) if self.hybband is not None else None
            hybbandpath, hybbandvasprunxml = checkpath(hybbandpathjoin, "vasprun.xml")

            absorbpathjion = os.path.join(self.workdir, material.strip(), self.optical) if self.optical is not None else None
            absorbpath, absorbvasprunxml = checkpath(absorbpathjion, "vasprun.xml")

            esppathjion = os.path.join(self.workdir, material.strip(), self.esp) if self.esp is not None else None
            esppath, locpotfile = checkpath(esppathjion, "LOCPOT")
            esppath, espoutcarfile = checkpath(esppathjion, "OUTCAR")

            if self.bandplot:
                if not os.path.exists(os.path.join(self.workdir, "BAND_dir")):
                    os.makedirs(os.path.join(self.workdir, "BAND_dir"))

            if self.pmg:
                self.info_dict[material] = self.pmggetinfo(material, stfile, scfvasprunxml, scfoutcarfile,
                                                           bandvasprunxml,
                                                           hybbandvasprunxml, absorbvasprunxml, locpotfile,
                                                           espoutcarfile)
            elif self.vaspkit:
                self.info_dict[material] = self.vaspkitinfo(material, stfile, scfoutcarfile, scfoszicarfile, bandpath,
                                                            hybbandpath, absorbpath, esppath)
            else:
                self.info_dict[material] = self.vaspkitinfo(material, stfile, scfoutcarfile, scfoszicarfile,
                                                            None, None, None, None)

        with open(self.infosavefile, "w") as f:
            json.dump(self.info_dict, f, cls=NumpyEncoder)

    def getInfoBi(self, monodict, bimondict=None, natomlayer1=3, hybrid=False, overwrite=True):
        """
        Organize the dictionary for bilayer materials.

        Args:
            monodict: dict, the dictionary of monolayer materials. Two types of dictionaries are supported:
                1.  {"1MoS2": {"label": SMoS, ...}}, the directory name as the key and contains "label" information.
                2.  {"SMoS": {"material": "1MoS2", ...}}, the label as the key and contains "material" information.
            bimondict: dict, the dictionary of monolayer materials contained in bilayer materials.
            natomlayer1: int, the number of atoms in the first layer.
            hybrid: bool, whether to plot the hybrid band structure.
            overwrite: bool, whether to overwrite the existing dictionary of bilayer materials.
        """
        if self.bandplot:
            if not os.path.exists(os.path.join(self.workdir, "LayerBand_dir")):
                os.mkdir(os.path.join(self.workdir, "LayerBand_dir"))

        monodicttype = None
        bimondictload = None
        if os.path.exists(os.path.join(self.workdir, monodict)):
            with open(os.path.join(self.workdir, monodict), "r") as f:
                mono_dict = json.load(f)
            monoinfokeys = list(mono_dict.values())[0].keys()
            if "label" in monoinfokeys:
                monodicttype = "label"
            elif "material" in monoinfokeys:
                monodicttype = "material"
        if monodicttype is None:
            warnings.warn(f"The monolayer info_dict file {monodict} does not exist.")
        else:
            if bimondict is not None and os.path.exists(os.path.join(self.workdir, bimondict)):
                with open(os.path.join(self.workdir, bimondict), "r") as f:
                    bimondictload = json.load(f)
                self.infomono_dict = mono_dict
            else:
                self.infomono_dict = {}
                if monodicttype == "material":
                    for mono, monoinfo in mono_dict.items():
                        if mono not in self.infomono_dict:
                            self.infomono_dict[mono] = monoinfo
                            if "elements" in monoinfo:
                                mono_inv = "".join(monoinfo["elements"][::-1])
                                if mono_inv != mono:
                                    self.infomono_dict[mono_inv] = monoinfo
                else:
                    for mono, monoinfo in mono_dict.items():
                        label = monoinfo["label"]
                        if label not in self.infomono_dict:
                            self.infomono_dict[label] = monoinfo
                            self.infomono_dict[label]["material"] = mono
                            if "elements" in monoinfo:
                                mono_inv = "".join(monoinfo["elements"][::-1])
                                if mono_inv != label:
                                    self.infomono_dict[mono_inv] = monoinfo

        infokeys = list(self.info_dict.values())[0].keys()
        required_keys = ["elements", "natoms", "atom_indices", "zdiff", "Efermi", "E0"]

        for key_ in required_keys:
            if key_ not in infokeys:
                raise ValueError(f"{key_} is not in the info_dict, please re-initialize the info_dict.")

        self.infobi_dict = {}

        def sumpro(material_data, materialid, natomlayer1, prefix=""):
            if "band" not in material_data or prefix + "vbmpro" not in material_data["band"]:
                self.infobi_dict[materialid][prefix + "vbmprosum"].append((None, None))
                self.infobi_dict[materialid][prefix + "cbmprosum"].append((None, None))
                return None, None, None, None

            vbmpro = material_data["band"][prefix + "vbmpro"]
            cbmpro = material_data["band"][prefix + "cbmpro"]
            vbmpro_l1 = sum(val_b for val_a, val_b in zip(material_data["atom_indices"], vbmpro) if val_a < natomlayer1)
            vbmpro_l2 = sum(
                val_b for val_a, val_b in zip(material_data["atom_indices"], vbmpro) if val_a >= natomlayer1)
            cbmpro_l1 = sum(val_b for val_a, val_b in zip(material_data["atom_indices"], cbmpro) if val_a < natomlayer1)
            cbmpro_l2 = sum(
                val_b for val_a, val_b in zip(material_data["atom_indices"], cbmpro) if val_a >= natomlayer1)
            self.infobi_dict[materialid][prefix + "vbmprosum"].append((vbmpro_l1, vbmpro_l2))
            self.infobi_dict[materialid][prefix + "cbmprosum"].append((cbmpro_l1, cbmpro_l2))
            return vbmpro_l1, vbmpro_l2, cbmpro_l1, cbmpro_l2

        for material, material_data in self.info_dict.items():
            l1E0, l2E0, ismetal1, ismetal2 = None, None, None, None
            vbmpro_l1, vbmpro_l2, cbmpro_l1, cbmpro_l2 = None, None, None, None
            direct_vbmpro_l1, direct_vbmpro_l2, direct_cbmpro_l1, direct_cbmpro_l2 = None, None, None, None
            ismetalb = None

            material_parts = material.strip().split("/")
            material_info = material_parts[0].split("_")

            materialid = material_info[0]
            mismatch = material_info[1] if len(material_info) > 1 else None
            natomlayer1 = int(material_info[2]) if len(material_info) > 2 else natomlayer1
            stackmode = "/".join(material_parts[1:]).rstrip("/") if len(material_parts) > 1 else None

            if materialid not in self.infobi_dict:
                self.infobi_dict[materialid] = {
                    "elements": material_data["elements"],
                    "natoms": material_data["natoms"],
                    "atom_indices": material_data["atom_indices"],
                    "stackmode": [stackmode],
                    "Efermi": [material_data["Efermi"]],
                    "E0": [material_data["E0"]]
                }

                if "mag" in material_data:
                    self.infobi_dict[materialid]["mag"] = [material_data["mag"]]

                if mismatch:
                    self.infobi_dict[materialid]["mismatch"] = [mismatch]
                # if natomlayer1:
                #     self.infobi_dict[materialid]["natomlayer1"] = [natomlayer1]

                if natomlayer1:
                    if bimondictload:
                        bl1label = bimondictload[material_parts[0]][0]
                        bl2label = bimondictload[material_parts[0]][1]
                    else:
                        bl1label = "".join(material_data["elements"][:natomlayer1])
                        bl2label = "".join(material_data["elements"][natomlayer1:])
                    self.infobi_dict[materialid]["natomlayer1"] = [natomlayer1]
                    self.infobi_dict[materialid]["blabel"] = material_data["label"]
                    self.infobi_dict[materialid]["bl1label"] = bl1label
                    self.infobi_dict[materialid]["bl2label"] = bl2label
                    self.infobi_dict[materialid]["dinterlayer"] = [material_data["zdiff"][-natomlayer1]]

                    if monodicttype is not None:
                        if bl1label in self.infomono_dict and bl2label in self.infomono_dict:
                            l1E0 = self.infomono_dict[bl1label]["E0"]
                            l2E0 = self.infomono_dict[bl2label]["E0"]
                            self.infobi_dict[materialid]["l1info"] = {"E0": l1E0}
                            self.infobi_dict[materialid]["l2info"] = {"E0": l2E0}
                            self.infobi_dict[materialid]["Ef"] = [material_data["E0"] - l1E0 - l2E0]
                            la1 = self.infomono_dict[bl1label]["abc"][0]
                            la2 = self.infomono_dict[bl2label]["abc"][0]
                            self.infobi_dict[materialid]["mismatchopt"] = 2 * abs(la1 - la2) / (la1 + la2)
                            if "band" in monoinfokeys:
                                bandgapl1 = self.infomono_dict[bl1label]["band"]["bandgap"]
                                bandgapl2 = self.infomono_dict[bl2label]["band"]["bandgap"]
                                ismetal1 = self.infomono_dict[bl1label]["band"]["is_metals"]
                                ismetal2 = self.infomono_dict[bl2label]["band"]["is_metals"]
                                self.infobi_dict[materialid]["l1info"].update({
                                    "bandgap": bandgapl1,
                                    "is_metals": ismetal1
                                })
                                self.infobi_dict[materialid]["l2info"].update({
                                    "bandgap": bandgapl2,
                                    "is_metals": ismetal2
                                })

                if "band" in infokeys:
                    ismetalb = material_data["band"]["is_metals"]
                    bandgap = material_data["band"]["bandgap"]
                    self.infobi_dict[materialid].update({
                        "bandgap": [bandgap],
                        "direct": [material_data["band"]["direct"]],
                        "is_metals": [ismetalb],
                        "vbmprosum": [],
                        "cbmprosum": [],
                        "bandalignment": [],
                        "deltabandgap": [],
                        "direct_vbmprosum": [],
                        "direct_cbmprosum": [],
                        "directbandalignment": []
                    })
                    vbmpro_l1, vbmpro_l2, cbmpro_l1, cbmpro_l2 = sumpro(material_data, materialid, natomlayer1)
                    if "directbandgap" in material_data["band"]:
                        if material_data["band"]["directbandgap"]:
                            directbandgap = material_data["band"]["directbandgap"]["1"]["value"]
                            if "-1" in material_data["band"]["directbandgap"]:
                                directbandgap = min(directbandgap,
                                                    material_data["band"]["directbandgap"]["-1"]["value"])
                            self.infobi_dict[materialid]["deltabandgap"] = [directbandgap - bandgap]
                            direct_vbmpro_l1, direct_vbmpro_l2, direct_cbmpro_l1, direct_cbmpro_l2 = sumpro(
                                material_data, materialid, natomlayer1, prefix="direct_")
                        else:
                            self.infobi_dict[materialid]["deltabandgap"] = [1.0]
                    else:
                        self.infobi_dict[materialid]["deltabandgap"] = [1.0]
            else:
                self.infobi_dict[materialid]["stackmode"].append(stackmode)
                self.infobi_dict[materialid]["Efermi"].append(material_data["Efermi"])
                self.infobi_dict[materialid]["E0"].append(material_data["E0"])
                self.infobi_dict[materialid]["dinterlayer"].append(material_data["zdiff"][-natomlayer1])
                if mismatch:
                    self.infobi_dict[materialid]["mismatch"].append(mismatch)
                if natomlayer1:
                    self.infobi_dict[materialid]["natomlayer1"].append(natomlayer1)
                if "mag" in material_data:
                    self.infobi_dict[materialid]["mag"].append(material_data["mag"])

                if "band" in infokeys:
                    ismetalb = material_data["band"]["is_metals"]
                    ismetal1 = self.infomono_dict[self.infobi_dict[materialid]["bl1label"]]["band"]["is_metals"]
                    ismetal2 = self.infomono_dict[self.infobi_dict[materialid]["bl2label"]]["band"]["is_metals"]
                    bandgap = material_data["band"]["bandgap"]
                    self.infobi_dict[materialid]["bandgap"].append(bandgap)
                    if "directbandgap" in material_data["band"]:
                        if material_data["band"]["directbandgap"]:
                            directbandgap = material_data["band"]["directbandgap"]["1"]["value"]
                            if "-1" in material_data["band"]["directbandgap"]:
                                directbandgap = min(directbandgap,
                                                    material_data["band"]["directbandgap"]["-1"]["value"])
                            self.infobi_dict[materialid]["deltabandgap"].append(directbandgap - bandgap)
                            direct_vbmpro_l1, direct_vbmpro_l2, direct_cbmpro_l1, direct_cbmpro_l2 = sumpro(
                                material_data, materialid, natomlayer1, prefix="direct_")
                        else:
                            self.infobi_dict[materialid]["deltabandgap"].append(1.0)
                    else:
                        self.infobi_dict[materialid]["deltabandgap"].append(1.0)
                    self.infobi_dict[materialid]["direct"].append(material_data["band"]["direct"])
                    self.infobi_dict[materialid]["is_metals"].append(ismetalb)
                    vbmpro_l1, vbmpro_l2, cbmpro_l1, cbmpro_l2 = sumpro(material_data, materialid, natomlayer1)

                if "Ef" in self.infobi_dict[materialid].keys():
                    l1E0 = self.infobi_dict[materialid]["l1info"]["E0"]
                    l2E0 = self.infobi_dict[materialid]["l2info"]["E0"]
                    self.infobi_dict[materialid]["Ef"].append(material_data["E0"] - l1E0 - l2E0)

            if "band" in infokeys:
                if vbmpro_l1 is not None:
                    bandalignment = self.bandalignment(ismetal1, ismetal2, ismetalb,
                                                       vbmpro_l1, vbmpro_l2, cbmpro_l1, cbmpro_l2)
                    self.infobi_dict[materialid]["bandalignment"].append(bandalignment)
                else:
                    self.infobi_dict[materialid]["bandalignment"].append(None)
                if "directbandalignment" in self.infobi_dict[materialid].keys():
                    if direct_vbmpro_l1 is not None:
                        directbandalignment = self.bandalignment(ismetal1, ismetal2, ismetalb,
                                                                 direct_vbmpro_l1, direct_vbmpro_l2, direct_cbmpro_l1,
                                                                 direct_cbmpro_l2)
                        self.infobi_dict[materialid]["directbandalignment"].append(directbandalignment)
                    else:
                        self.infobi_dict[materialid]["directbandalignment"].append(None)

            if self.bandplot and natomlayer1 and material_data["atom_indices"]:
                bandpath = os.path.join(self.workdir, material.strip(), self.band)
                plotBS(bandpath=bandpath, savename=material.strip(),
                       savepath=os.path.join(self.workdir, "LayerBand_dir"),
                       indices=material_data["atom_indices"], natomlayer1=natomlayer1,
                       erange=self.banderange, elempro=self.elempro, hybrid=hybrid, imag_format=self.imag_format,
                       pmg=self.pmg, vaspkit=self.vaspkit, BSDOSPlotter=self.BSDOSPlotter,
                       BSPlotter=self.BSPlotter, Vasprun=self.Vasprun).plotbsdos()

        bisavename = "biinfo"
        bisavepath = os.path.join(self.workdir, f"{bisavename}.json")
        bisavenu = 0
        if not overwrite:
            while os.path.exists(bisavepath):
                bisavenu += 1
                bisavename = f"biinfo{bisavenu}"
                bisavepath = os.path.join(self.workdir, f"{bisavename}.json")
        with open(bisavepath, "w") as f:
            json.dump(self.infobi_dict, f, cls=NumpyEncoder)

        bisimple_dict = {}
        infokeys = self.infobi_dict[list(self.infobi_dict.keys())[0]].keys()
        for materialid, material_data in self.infobi_dict.items():
            bisimple_dict[materialid] = {}
            # if "mismatch" in infokeys:
            #     bisimple_dict[materialid]["mismatch"] = material_data["mismatch"]
            if "mismatchopt" in infokeys:
                bisimple_dict[materialid]["mismatchopt"] = material_data["mismatchopt"]
            bisimple_dict[materialid].update({
                "stackmode": material_data["stackmode"],
                "Ef": material_data["Ef"],
                "dinterlayer": material_data["dinterlayer"]
            })
            if "bandalignment" in infokeys:
                bisimple_dict[materialid].update({
                    "bandgap": material_data["bandgap"],
                    "direct": material_data["direct"],
                    "bandalignment": material_data["bandalignment"]
                })
            if "directbandalignment" in infokeys:
                bisimple_dict[materialid]["directbandalignment"] = material_data["directbandalignment"]
            if "deltabandgap" in infokeys:
                bisimple_dict[materialid]["deltabandgap"] = material_data["deltabandgap"]
        with open(os.path.join(self.workdir, f"{bisavename}-simple.json"), "w") as f:
            json.dump(bisimple_dict, f, indent=4, cls=NumpyEncoder)

    def bandalignment(self, ismetal1, ismetal2, ismetalb, vbmpro_l1, vbmpro_l2, cbmpro_l1, cbmpro_l2):
        """
        Identifying an energy band alignment formed by the stacking of two semiconductors.
        0: Metal; 1: Type I; 2: Type II; 3: Type III.

        Args:
            l1ismetal (bool): Whether the first layer is metal.
            l2ismetal (bool): Whether the second layer is metal.
            bismetal (bool): Whether the bilayer is metal.
            vbmpro (float): Contributions of the atoms of the different layers in the VBM.
            cbmpro (float): Contributions of the atoms of the different layers in the CBM.
        """
        if ismetalb:
            if ismetal1 or ismetal2:
                return None
            if (vbmpro_l1 >= vbmpro_l2 and cbmpro_l1 >= cbmpro_l2) or (
                    vbmpro_l1 <= vbmpro_l2 and cbmpro_l1 <= cbmpro_l2):
                return 0
            else:
                return 3
        if (vbmpro_l1 >= vbmpro_l2 and cbmpro_l1 >= cbmpro_l2) or (vbmpro_l1 <= vbmpro_l2 and cbmpro_l1 <= cbmpro_l2):
            return 1
        else:
            return 2

    def pmggetinfo(self, material, stfile, scfvasprunxml, scfoutcarfile, bandvasprunxml, hybbandvasprunxml,
                   absorbvasprunxml, locpotfile, espoutcarfile):
        propertiesinfo = {}
        st = self.Structure.from_file(stfile)
        try:
            elements = st.labels
        except AttributeError:
            elements = [site.specie.symbol for site in st.sites]
        lattice = st.lattice.matrix
        abc = st.lattice.abc
        cart_coords = st.cart_coords
        indices = np.argsort(st.cart_coords[:, 2])
        ranking = np.empty_like(indices)
        ranking[indices] = np.arange(len(cart_coords))
        cart_coords = cart_coords[indices]
        elements = [elements[i] for i in indices]
        natoms = len(elements)
        zdiff = np.diff(cart_coords[:, 2])
        try:
            spacegroup = st.get_space_group_info()
        except AttributeError:
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            spg_analyzer = SpacegroupAnalyzer(st)
            spacegroup_data = spg_analyzer.get_symmetry_dataset()
            spacegroup_symbol = spacegroup_data['international'] if 'international' in spacegroup_data else None
            spacegroup_number = spacegroup_data['number'] if 'number' in spacegroup_data else None
            spacegroup = (spacegroup_symbol, spacegroup_number)
        propertiesinfo.update({
            "elements": elements,
            "label": "".join(elements),
            "natoms": natoms,
            "atom_indices": ranking,
            "lattice": lattice.tolist(),
            "abc": abc,
            "cart_coords": cart_coords.tolist(),
            "zdiff": zdiff.tolist(),
            "spacegroup": spacegroup
        })

        if scfvasprunxml:
            vasprun = self.Vasprun(scfvasprunxml)
            propertiesinfo["Efermi"] = vasprun.eigenvalue_band_properties[1]
            propertiesinfo["E0"] = vasprun.final_energy
            if self.mag and scfoutcarfile:
                outcar = self.Outcar(scfoutcarfile)
                propertiesinfo["Magnetization_x"] = outcar.magnetization
                scfoszicarfile = os.path.join(os.path.dirname(scfoutcarfile), "OSZICAR")
                mag_ = run_subprocess(f"tail -1 {scfoszicarfile} | awk '{{print $10}}'")
                propertiesinfo["mag"] = float(mag_) if mag_ is not None else None

        def getbandinfo(material, vasprunxml, hybrid=False):
            bandbasicinfo = {}
            force_hybrid_mode = True if hybrid else False
            bs = self.Vasprun(vasprunxml, parse_potcar_file=False).get_band_structure(line_mode=True,
                                                                                      force_hybrid_mode=force_hybrid_mode)
            if not bs.is_metal():
                vbm_info = bs.get_vbm()
                cbm_info = bs.get_cbm()
                vbm_info_dict, cbm_info_dict = {}, {}
                vbm_spin, vbm_band_index, cbm_spin, cbm_band_index = None, None, None, None
                for spin, band_indices in vbm_info['band_index'].items():
                    vbm_spin, vbm_band_index = str(spin), band_indices
                for spin, band_indices in cbm_info['band_index'].items():
                    cbm_spin, cbm_band_index = str(spin), band_indices
                vbm_info_dict["spin"], cbm_info_dict["spin"] = vbm_spin, cbm_spin
                vbm_info_dict["band_index"], cbm_info_dict["band_index"] = vbm_band_index, cbm_band_index
                vbm_info_dict["kpoint_index"], cbm_info_dict["kpoint_index"] = vbm_info['kpoint_index'], cbm_info[
                    'kpoint_index']

                bandbasicinfo.update({
                    "bandgap": bs.get_band_gap()['energy'],
                    "direct": bs.get_band_gap()['direct'],
                    "transition": bs.get_band_gap()['transition'],
                    "is_metals": bs.is_metal(),
                    "vbm": vbm_info['energy'],
                    "cbm": cbm_info['energy'],
                    "vbm_index": vbm_info_dict,
                    "cbm_index": cbm_info_dict,
                    "directbandgap": {str(key): value for key, value in bs.get_direct_band_gap_dict().items()}
                })
            else:
                bandbasicinfo.update({
                    "bandgap": 0.0,
                    "direct": False,
                    "transition": None,
                    "is_metals": bs.is_metal(),
                    "vbm": None,
                    "cbm": None,
                    "vbm_index": None,
                    "cbm_index": None,
                    "directbandgap": None
                })
            if self.nlayer == 2:
                vbm_index = bandbasicinfo["vbm_index"]
                cbm_index = bandbasicinfo["cbm_index"]
                procarpath = os.path.dirname(vasprunxml)
                if vbm_index is not None and cbm_index is not None:
                    vbmpro, cbmpro = readPROCAR(procarpath, vbm_index, cbm_index, self.weightsave)
                    bandbasicinfo.update({"vbmpro": vbmpro, "cbmpro": cbmpro})
                if bandbasicinfo["directbandgap"] is not None:
                    direct_vbm_info_dict, direct_cbm_info_dict = {}, {}
                    direct_spin = "1"
                    if "-1" in bandbasicinfo["directbandgap"]:
                        direct_spin = "-1" if bandbasicinfo["directbandgap"]["-1"]["value"] < \
                                              bandbasicinfo["directbandgap"]["1"]["value"] else "1"
                    direct_vbm_info_dict["spin"], direct_cbm_info_dict["spin"] = direct_spin, direct_spin
                    direct_vbm_info_dict["band_index"], direct_cbm_info_dict["band_index"] = \
                        bandbasicinfo["directbandgap"][direct_spin]["band_indices"][0], \
                        bandbasicinfo["directbandgap"][direct_spin]["band_indices"][1]
                    direct_vbm_info_dict["kpoint_index"] = \
                        direct_cbm_info_dict["kpoint_index"] = bandbasicinfo["directbandgap"][direct_spin][
                        "kpoint_index"]
                    direct_vbmpro, direct_cbmpro = readPROCAR(procarpath, direct_vbm_info_dict, direct_cbm_info_dict,
                                                              self.weightsave)
                    bandbasicinfo.update({"direct_vbmpro": direct_vbmpro, "direct_cbmpro": direct_cbmpro})

            if self.bandplot:
                plotBS(bandpath=os.path.dirname(vasprunxml), savename=material.strip(),
                       savepath=os.path.join(self.workdir, "BAND_dir"),
                       erange=self.banderange, elempro=self.elempro, hybrid=hybrid, imag_format=self.imag_format,
                       pmg=self.pmg, vaspkit=self.vaspkit, BSDOSPlotter=self.BSDOSPlotter,
                       BSPlotter=self.BSPlotter, Vasprun=self.Vasprun).plotbsdos()
            return bandbasicinfo

        bandinfo = getbandinfo(material=material, vasprunxml=bandvasprunxml) if bandvasprunxml is not None else None
        hybbandinfo = getbandinfo(material=material, vasprunxml=hybbandvasprunxml,
                                  hybrid=True) if hybbandvasprunxml is not None else None

        absorbinfo = None
        if absorbvasprunxml is not None:
            absorbinfo = {}
            absorbinfo.update({"dielectric": self.Vasprun(absorbvasprunxml).dielectric})

        espinfo = None
        if locpotfile is not None and espoutcarfile is not None:
            espinfo = {}
            locpot = self.Locpot.from_file(locpotfile)
            z_potential = locpot.get_average_along_axis(2)
            vp1 = np.mean(z_potential[-int(0.05 * len(z_potential)):])
            vp2 = np.mean(z_potential[:int(0.05 * len(z_potential))])
            vacuum_potential = vp1 if vp1 > vp2 else vp2
            outcar = self.Outcar(espoutcarfile)
            fermi_level = outcar.efermi
            work_function = vacuum_potential - fermi_level

            espinfo.update({
                "coords_z": locpot.get_axis_grid(2),
                "potential_z": z_potential,
                "vacuum_potential": vacuum_potential,
                "work_function": work_function,
                "vacuum_potential_dict": {"vp_left": vp2, "vp_right": vp1, "vp_delta": vp2 - vp1}
            })

        infokey = ["band", "hybband", "absorption", "esp"]
        for i_, info_ in enumerate([bandinfo, hybbandinfo, absorbinfo, espinfo]):
            if info_:
                propertiesinfo.update({infokey[i_]: info_})
        return propertiesinfo

    def vaspkitinfo(self, material, stfile, scfoutcarfile, scfoszicarfile, bandpath, hybbandpath, absorbpath, esppath):
        propertiesinfo = {}
        st = StructureFile(stfile)
        elements = st.elements
        natoms = len(elements)
        lattice = st.lattice
        abc = st.abc
        cart_coords = st.positions
        zdiff = st.zdiff
        indices = st.indices
        ranking = np.empty_like(indices)
        ranking[indices] = np.arange(len(cart_coords))
        spacegroup = None
        if self.vaspkit:
            matpath = os.path.dirname(stfile)
            if os.path.exists(os.path.join(matpath, "symmetry.log")):
                os.remove(os.path.join(matpath, "symmetry.log"))
            os.system(f"(cd {matpath} && vaspkit -task 601 >> symmetry.log)")
            sg = run_subprocess(f"grep 'Space Group' {matpath}/symmetry.log" + " | tail -n 1 | awk '{print $3}'")
            isg = run_subprocess(f"grep 'International' {matpath}/symmetry.log" + " | tail -n 1 | awk '{print $2}'")
            spacegroup = (isg, sg)
        propertiesinfo.update({
            "elements": elements,
            "label": "".join(elements),
            "natoms": natoms,
            "atom_indices": ranking,
            "lattice": lattice.tolist(),
            "abc": abc,
            "cart_coords": cart_coords.tolist(),
            "zdiff": zdiff.tolist()
        })
        if spacegroup:
            propertiesinfo.update({"spacegroup": spacegroup})

        if scfoutcarfile:
            propertiesinfo["Efermi"] = float(
                run_subprocess(f"grep 'E-fermi' {scfoutcarfile}" + " | tail -n 1 | awk '{print $3}'"))
        if scfoszicarfile:
            propertiesinfo["E0"] = float(run_subprocess(f"tail -1 {scfoszicarfile}" + " | awk '{print $5}'"))
        if self.mag:
            mag_ = run_subprocess(f"tail -1 {scfoszicarfile} | awk '{{print $10}}'")
            propertiesinfo["mag"] = float(mag_) if mag_ is not None else None

        def getbandinfo(material, path, hybrid=False, dim=2):
            bandbasicinfo = {}
            if hybrid:
                if not os.path.exists(os.path.join(path, "KPATH.in")):
                    assert dim in [1, 2, 3], "The dimension of structure should be 1, 2 or 3."
                    os.system(f"(cd {path} && vaspkit -task 30{dim} >/dev/null 2>&1)")
                os.system(f"(cd {path} && vaspkit -task 252 >/dev/null 2>&1)")
            else:
                os.system(f"(cd {path} && vaspkit -task 211 >/dev/null 2>&1)")

            if not os.path.exists(os.path.join(path, "BAND_GAP")):
                return None

            spin = run_subprocess(f"grep 'Spin' {path}/BAND_GAP | awk '{{print $3}}'")
            bandgap_up = float(run_subprocess(f"grep 'Gap' {path}/BAND_GAP | awk '{{print $4}}'"))
            vbm_up = float(run_subprocess(f"grep 'Eigenvalue of VBM' {path}/BAND_GAP | awk '{{print $5}}'"))
            cbm_up = float(run_subprocess(f"grep 'Eigenvalue of CBM' {path}/BAND_GAP | awk '{{print $5}}'"))
            vbm_info_dict, cbm_info_dict = {}, {}

            if spin:
                spinbandgap = {}
                bandgap_down = float(run_subprocess(f"grep 'Gap' {path}/BAND_GAP | awk '{{print $5}}'"))
                vbm_down = float(run_subprocess(f"grep 'Eigenvalue of VBM' {path}/BAND_GAP | awk '{{print $6}}'"))
                cbm_down = float(run_subprocess(f"grep 'Eigenvalue of CBM' {path}/BAND_GAP | awk '{{print $6}}'"))

                suffixlist = [" (UP)", " (DW)"]
                for sufi, suffix in enumerate(["1", "-1"]):
                    spinbandgap[suffix] = {
                        "bandgap": bandgap_up if suffix == "1" else bandgap_down,
                        "vbm": vbm_up if suffix == "1" else vbm_down,
                        "cbm": cbm_up if suffix == "1" else cbm_down
                    }
                    vbmlocation, cbmlocation = get_vbm_cbm_locations(path, suffixlist[sufi])
                    spinbandgap[suffix]["direct"] = is_direct_band_gap(vbmlocation, cbmlocation)

                vbmlocation, cbmlocation = get_vbm_cbm_locations(path)
                vbm_band_index, cbm_band_index, vbmspin, cbmspin = get_bm_spin(path)

                bandgap = float(run_subprocess(f"grep 'Gap' {path}/BAND_GAP" + " | awk '{print $6}'"))
                bandbasicinfo["bandgap"] = bandgap
                bandbasicinfo["direct"] = is_direct_band_gap(vbmlocation, cbmlocation)
                bandbasicinfo["is_metals"] = bandgap <= 0
                bandbasicinfo["vbm"] = float(
                    run_subprocess(f"grep 'Eigenvalue of VBM' {path}/BAND_GAP | awk '{{print $7}}'"))
                bandbasicinfo["cbm"] = float(
                    run_subprocess(f"grep 'Eigenvalue of CBM' {path}/BAND_GAP | awk '{{print $7}}'"))
                bandbasicinfo["vbm_index"] = {"spin": vbmspin, "band_index": vbm_band_index,
                                              "kpoint_coord": vbmlocation}
                bandbasicinfo["cbm_index"] = {"spin": cbmspin, "band_index": cbm_band_index,
                                              "kpoint_coord": cbmlocation}
                bandbasicinfo["spinbandgap"] = spinbandgap
            else:
                vbm_band_index = [
                    int(run_subprocess(f"grep 'HOMO & LUMO Bands' {path}/BAND_GAP | awk '{{print $5}}'")) - 1]
                cbm_band_index = [
                    int(run_subprocess(f"grep 'HOMO & LUMO Bands' {path}/BAND_GAP | awk '{{print $6}}'")) - 1]
                vbmlocation, cbmlocation = get_vbm_cbm_locations(path, "")
                vbm_info_dict.update({"spin": "1", "band_index": vbm_band_index, "kpoint_coord": vbmlocation})
                cbm_info_dict.update({"spin": "1", "band_index": cbm_band_index, "kpoint_coord": cbmlocation})

                bandbasicinfo.update({
                    "bandgap": bandgap_up,
                    "direct": run_subprocess(f"grep 'Character' {path}/BAND_GAP | awk '{{print $3}}'") == "Direct",
                    "is_metals": bandgap_up <= 0,
                    "vbm": vbm_up,
                    "cbm": cbm_up,
                    "vbm_index": vbm_info_dict,
                    "cbm_index": cbm_info_dict
                })

            if self.directbandgap and self.Vasprun is not None:
                force_hybrid_mode = True if hybrid else False
                vasprunxml = os.path.join(path, "vasprun.xml")
                bs = self.Vasprun(vasprunxml, parse_potcar_file=False).get_band_structure(line_mode=True,
                                                                                          force_hybrid_mode=force_hybrid_mode)
                if not bs.is_metal():
                    bandbasicinfo.update({"directbandgap": {str(key): value for key, value in
                                                            bs.get_direct_band_gap_dict().items()}})

            if self.nlayer == 2:
                vbm_index = bandbasicinfo["vbm_index"]
                cbm_index = bandbasicinfo["cbm_index"]
                if vbm_index is not None and cbm_index is not None:
                    vbmpro, cbmpro = readPROCAR(path, vbm_index, cbm_index, self.weightsave)
                    bandbasicinfo.update({"vbmpro": vbmpro, "cbmpro": cbmpro})
                if "directbandgap" in bandbasicinfo and bandbasicinfo["directbandgap"] is not None:
                    direct_vbm_info_dict, direct_cbm_info_dict = {}, {}
                    direct_spin = "1"
                    if "-1" in bandbasicinfo["directbandgap"]:
                        direct_spin = "-1" if bandbasicinfo["directbandgap"]["-1"]["value"] < \
                                              bandbasicinfo["directbandgap"]["1"]["value"] else "1"
                    direct_vbm_info_dict["spin"], direct_cbm_info_dict["spin"] = direct_spin, direct_spin
                    direct_vbm_info_dict["band_index"], direct_cbm_info_dict["band_index"] = \
                        bandbasicinfo["directbandgap"][direct_spin]["band_indices"][0], \
                        bandbasicinfo["directbandgap"][direct_spin]["band_indices"][1]
                    direct_vbm_info_dict["kpoint_index"] = \
                        direct_cbm_info_dict["kpoint_index"] = bandbasicinfo["directbandgap"][direct_spin][
                        "kpoint_index"]
                    direct_vbmpro, direct_cbmpro = readPROCAR(path, direct_vbm_info_dict,
                                                              direct_cbm_info_dict, self.weightsave)
                    bandbasicinfo.update({"direct_vbmpro": direct_vbmpro, "direct_cbmpro": direct_cbmpro})

            if self.bandplot:
                plotBS(bandpath=path, savename=material.strip(), savepath=os.path.join(self.workdir, "BAND_dir"),
                       erange=self.banderange, elempro=self.elempro, hybrid=hybrid, imag_format=self.imag_format,
                       pmg=self.pmg, vaspkit=self.vaspkit, BSDOSPlotter=self.BSDOSPlotter,
                       BSPlotter=self.BSPlotter, Vasprun=self.Vasprun).plotbsdos()
            return bandbasicinfo

        bandinfo = getbandinfo(material=material, path=bandpath, dim=self.dim) if bandpath is not None else None
        hybbandinfo = getbandinfo(material=material, path=hybbandpath, hybrid=True,
                                  dim=self.dim) if hybbandpath is not None else None

        absorbinfo = None
        if absorbpath is not None:
            absorbinfo = {}
            os.system(f"echo -e '71{self.dim - 2}\\n{self.energyunit}' | (cd {absorbpath} && vaspkit >/dev/null 2>&1)")

            if os.path.exists(os.path.join(absorbpath, "REAL.in")):
                real_data = np.loadtxt(os.path.join(absorbpath, "REAL.in"))
                imag_data = np.loadtxt(os.path.join(absorbpath, "IMAG.in"))
                frequencies = real_data[:, 0]
                dielectric_real = real_data[:, 1:]
                dielectric_imag = imag_data[:, 1:]
                dielectric = [frequencies.tolist(), dielectric_real.tolist(), dielectric_imag.tolist()]
                absorbinfo.update({"dielectric": dielectric})

            if os.path.exists(os.path.join(absorbpath, "ABSORPTION_2D.dat")):
                datapath = os.path.join(absorbpath, "ABSORPTION_2D.dat")
            elif os.path.exists(os.path.join(absorbpath, "ABSORPTION.dat")):
                datapath = os.path.join(absorbpath, "ABSORPTION.dat")
            else:
                datapath = None
            if datapath:
                absorption_data = np.loadtxt(datapath)
                energies = absorption_data[:, 0]
                absorption_coefficients = absorption_data[:, 1:]
                absorbinfo.update({
                    "energies": energies.tolist(),
                    "absorption_coefficients": absorption_coefficients.tolist()
                })

        espinfo = None
        if esppath is not None:
            espinfo = {}
            os.system(f"echo -e '426\\n3' | (cd {esppath} && vaspkit >/dev/null 2>&1)")
            if os.path.exists(os.path.join(esppath, "PLANAR_AVERAGE.dat")):
                esp_data = np.loadtxt(os.path.join(esppath, "PLANAR_AVERAGE.dat"))
                espinfo.update({
                    "coords_z": esp_data[:, 0].tolist(),
                    "potential_z": esp_data[:, 1].tolist()
                })
                vp1 = np.mean(esp_data[:, 1][-int(0.05 * len(esp_data[:, 1])):])
                vp2 = np.mean(esp_data[:, 1][:int(0.05 * len(esp_data[:, 1]))])
                os.system(f"echo -e '927' | (cd {esppath} && vaspkit > Vacuum.txt)")
                if os.path.exists(os.path.join(esppath, "Vacuum.txt")):
                    vacuum_potential = float(
                        run_subprocess(f"grep 'Vacuum Level' {esppath}/Vacuum.txt | tail -n 1 | awk '{{print $4}}"))
                    work_function = float(
                        run_subprocess(f"grep 'Work Function' {esppath}/Vacuum.txt | tail -n 1 | awk '{{print $4}}"))
                    espinfo.update({
                        "vacuum_potential": vacuum_potential,
                        "work_function": work_function,
                        "vacuum_potential_dict": {"vp_left": vp2, "vp_right": vp1, "vp_delta": vp2 - vp1}
                    })
            else:
                print(f"Cannot find the ESP.dat file in {esppath}.")
        infokey = ["band", "hybband", "absorption", "esp"]
        for i_, info_ in enumerate([bandinfo, hybbandinfo, absorbinfo, espinfo]):
            if info_:
                propertiesinfo.update({infokey[i_]: info_})
        return propertiesinfo

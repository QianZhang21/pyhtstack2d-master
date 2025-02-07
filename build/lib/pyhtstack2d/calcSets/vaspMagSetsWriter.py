import os
import json
# import shutil

import subprocess
import warnings

import numpy as np

from .Pu import INCARPu

subdict = None
SUB_PATH = os.path.join(os.path.expanduser("~"), ".config", ".PyHTStack2D.json")
if os.path.exists(SUB_PATH):
    with open(SUB_PATH, 'r') as f:
        subdict = json.load(f)

maginitpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'maginit.json')
with open(maginitpath, "r") as f:
    maginit = json.load(f)


def determine_mag_state(magnetizations, maglowerl):
    # Check if all magnetization values are close to zero
    if all(abs(m) < maglowerl for m in magnetizations):  # 1e-3 is a small threshold you can adjust
        if abs(sum(magnetizations)) <= maglowerl:
            return "NM"
        else:
            positive = any(m > (maglowerl * 0.8) for m in magnetizations)
            negative = any(m < -(maglowerl * 0.8) for m in magnetizations)
            if positive and negative:
                return "AFM"
            else:
                return "FM"

    # Check for signs of magnetizations
    positive = any(m > maglowerl for m in magnetizations)
    negative = any(m < -maglowerl for m in magnetizations)

    if positive and negative:
        return "AFM"  # Antiferromagnetic
    else:
        return "FM"  # Ferromagnetic


class MonoMagState:
    def __init__(self, maglist="maglist.txt", workdir=None, magpath=None, magnetfile="magnetizationx.txt", pmg=False,
                 maglowerl=0.5, nomagatom=["Cl", "F", "I", "Br"], natomlayer1=3, skipmagnetfile=False, pU=True,
                 magmom_setting=None, magatomlist=None, onlymagpu=False, mkfm=True, **kwargs):
        """
        Generation of input files for magnetic ground state determination.

        maglist: str
            The file path containing the path where the magnetic system is located, usually in the form of TXT. However, if the file is in the form of JSON (e.g., magdict.json), the file will be loaded as a dictionary to initialise the magdict.
            This file can be generated as 'maglist.txt' using the separateMag function of the GetMagEntropy class of the genInput module.
            Once it has been determined that all the materials are magnetic in the current directory, a 'maglist.txt' file can be created using a command such as the following:
            ---------------V1------------------
            #!/bin/bash -x

            rm -f maglist.txt
            multilevel=2
            path=$(pwd -P)
            find "$path" -mindepth $multilevel -type d | while read -r dir; do
                echo "$dir" >> maglist.txt
            done
            -----------------------------------
            ---------------V2------------------
            #!/bin/bash -x

            rm -f maglist.txt
            path=$(pwd -P)
            for file in `ls $path`; do
            if [ -d "$file" ]; then
            for file2 in `ls $file`; do
            materialpath="$file/$file2"
            if [ -d $materialpath ]; then
            echo $materialpath >> maglist.txt
            fi
            done
            fi
            done
            -----------------------------------
        magpath: str
            The path of the VASP calculation result to get the magnetic moments.
            magpath="scf/" means that the magnetic moments are in the OSZICAR of "scf" directory.
            magpath="-scf" means that the magnetic moments are in the OSZICAR-scf of material directory.
            If None, it means that the magnetic moments are in the OSZICAR of material directory.
        magnetfile: str
            The file path containing the magnetic moments.
            The parameter "LORBIT = 11" should be set in the INCAR file.
        pmg: bool
            Whether to use the pymatgen module to get the magnetic moments, when the magnetfile is not available.
        maglowerl: float
            Identifying magnetic atoms, the lower limit of the magnetic moment.
        nomagatom: list
            The list of atoms that are not magnetic.
        natomlayer1: int
            The number of atoms of the first layer.
        skipmagnetfile: bool
            For bilayer systems, whether to skip reading the magnetic moments and judge the magnetic atoms directly based on the monolayer information provided.
        pU: bool
            Whether to add the U parameter to the INCAR file.
        magmom_setting: dict
            The dictionary of the magnetic moments of the atoms.
        magatomlist: list
            When the skipmagnetfile is True, the list of magnetic atoms can be provided directly.
        onlymagpu: bool
            Whether to add the U parameter for magnetic atoms only.
        mkfm: bool
            Forcing to create FM folders.
        """
        self.nomagatom = nomagatom
        # self.maglist = maglist
        self.workdir = workdir
        if self.workdir is None:
            # self.workdir = os.getcwd()
            self.workdir = ""
        self.maglist = os.path.join(self.workdir, maglist)
        self.pU = pU
        self.magmom_setting = magmom_setting
        self.magatomlist = magatomlist
        self.onlymagpu = onlymagpu
        self.mkfm = mkfm

        self.kpointspath = "KPOINTS"
        self.magpath = magpath
        if self.magpath is None:
            self.outpath = "OUTCAR"
            self.pospath = "CONTCAR"
        elif "/" in self.magpath:
            self.outpath = magpath + "OUTCAR"
            self.pospath = magpath + "CONTCAR"
            self.kpointspath = magpath + "KPOINTS"
        elif "-" in self.magpath:
            self.outpath = f"OUTCAR{magpath}"
            self.pospath = f"CONTCAR{magpath}"
        else:
            raise ValueError("magpath should be None, or a string with '/' or '-', such as 'scf/' or '-scf'.")

        self.magnetfile = magnetfile
        self.pmg = pmg
        self.Outcar = None
        self.Structure = None
        if self.pmg:
            self.Outcar = __import__('pymatgen.io.vasp.outputs', fromlist=['Outcar']).Outcar
            self.Structure = __import__('pymatgen.core.structure', fromlist=['Structure']).Structure

        self.natomlayer1 = natomlayer1
        self.skipmagnetfile = skipmagnetfile
        self.maglowerl = maglowerl
        if self.maglist.endswith(".json"):
            with open(self.maglist, "r") as f:
                self.magdict = json.load(f)
        else:
            self.magdict = self.analyseMagnet()

        self.iniMonoMagstates(natomlayer1=self.natomlayer1)

        for keys_ in self.magdict.keys():
            keys_ = keys_.split("/")
            self.multilevel = len(keys_)
            break

        self.runsh = os.path.join(self.workdir, "run.sh")
        self.bicalmaterial = []

        self.kwargs = kwargs
        self.puset = {"u_setting": None, "openmixparam": True, "mixparam": None}
        for key in self.puset.keys():
            if key in self.kwargs:
                self.puset[key] = self.kwargs

    def iniMonoMagstates(self, natomlayer1=3):
        pass

    def getMagnet(self):
        """
        Get the magnetic moments.
        """
        grepfile = os.path.join(self.workdir, "grepMagx.sh")
        with open(grepfile, "w", newline="\n") as f:
            f.write("#!/bin/bash -x\n\n")
            f.write("rm -f magnetizationx.txt\n")
            f.write(f"for dir in `cat {self.maglist}`; do\n")
            f.write(f"    posfile=$dir/{self.pospath}\n")
            f.write("    element=$(cat $posfile | sed -n '6p' | sed 's/[ ][ ]*/-/g')\n")
            f.write("    elementnum=$(cat $posfile | sed -n '7p' | sed 's/[ ][ ]*/-/g')\n\n")
            f.write(f"    outfile=$dir/{self.outpath}\n")
            f.write("    startnu=$(grep -n '# of ion' $outfile | tail -n 1 | cut -d':' -f1)\n")
            f.write("    endnu=$(grep -n 'tot ' $outfile | tail -n 1 | cut -d':' -f1)\n")
            f.write("    nadd=2\n    nsub=1\n")
            f.write("    nbegin=$((startnu + nadd))\n")
            f.write("    nend=$((endnu - nsub))\n\n")
            f.write("    echo \"------- $dir $element $elementnum\" | sed 's/\\r//g' >> magnetizationx.txt\n")
            f.write("    for ((i=nbegin; i<nend; i++)); do\n")
            f.write("        mag=$(sed -n \"${i}p\" $outfile)\n")
            f.write("        echo $mag >> magnetizationx.txt\n")
            f.write("    done\n")
            f.write("    totmag=$(sed -n \"${endnu}p\" $outfile)\n")
            f.write("    echo $totmag >> magnetizationx.txt\n")
            f.write("done\n")
        os.system(f"bash {grepfile}")

    def analyseMagnet(self):
        """
        Analyse the magnetic moments.
        """
        magnetfile_ = os.path.join(self.workdir, self.magnetfile)
        if not self.pmg and not self.skipmagnetfile:
            if not os.path.exists(magnetfile_):
                self.getMagnet()
            try:
                with open(magnetfile_, "r") as f:
                    lines = f.readlines()
                magnet_dict = {}
                material = None
                atommag = []
                element = []
                for line in lines:
                    line = line.strip()
                    if "-------" in line:
                        if material:
                            magatom = []
                            for i, key_ in enumerate(magnet_dict[material].keys()):
                                if key_ != "tot":
                                    if abs(magnet_dict[material][key_][-1]) > self.maglowerl and \
                                            element[i] not in self.nomagatom:
                                        magatom.append(element[i])
                                    else:
                                        magatom.append(None)

                            magnet_dict[material]["magatom"] = magatom
                            magnet_dict[material]["magvalues"] = atommag[:-1]
                            magnet_dict[material]["element"] = element_ori
                        material = line.split()[1].strip()
                        element = []
                        element_ori = [x for x in line.split()[2].strip().split('-') if x]
                        elementnum = [x for x in line.split()[3].strip().split('-') if x]
                        for i, elem in enumerate(element_ori):
                            for j in range(int(elementnum[i])):
                                element.append(elem)
                        magnet_dict[material] = {}
                        atommag = []
                    else:
                        mag_values = line.split()
                        magnet_dict[material][mag_values[0]] = [float(mag) for mag in mag_values[1:]]
                        atommag.append(float(mag_values[-1]))
                magatom = []
                for i, key_ in enumerate(magnet_dict[material].keys()):
                    if key_ != "tot":
                        if abs(magnet_dict[material][key_][-1]) > self.maglowerl and \
                                element[i] not in self.nomagatom:
                            magatom.append(element[i])
                        else:
                            magatom.append(None)
                magnet_dict[material]["magatom"] = magatom
                magnet_dict[material]["element"] = element_ori
            except:
                raise ValueError("The magnetic moments file does not exist.")

        elif self.Structure and self.Outcar and not self.skipmagnetfile:
            magnet_dict = {}
            assert os.path.exists(self.maglist), "The magnetic list file does not exist."
            with open(self.maglist, "r") as f:
                lines = f.readlines()

            for line in lines:
                material = line.strip()
                outcar = os.path.join(self.workdir, material, self.outpath)
                poscar = os.path.join(self.workdir, material, self.pospath)
                if not os.path.exists(outcar) or not os.path.exists(poscar):
                    print(f"Warning: {outcar} or {poscar} does not exist.")
                    continue
                element_ori = np.loadtxt(poscar, skiprows=5, max_rows=1, dtype=str)
                magnet_dict[material] = {}
                st = self.Structure.from_file(poscar)
                try:
                    elements = st.labels
                except AttributeError:
                    elements = [site.specie.symbol for site in st.sites]

                magmom = self.Outcar(outcar).magnetization
                atommag = [magtot['tot'] for magtot in magmom]
                magatom = []
                for i, element in enumerate(elements):
                    if abs(atommag[i]) > self.maglowerl and element not in self.nomagatom:
                        magatom.append(element)
                    else:
                        magatom.append(None)
                magnet_dict[material]["magatom"] = magatom
                magnet_dict[material]["magvalues"] = atommag
                magnet_dict[material]["element"] = element_ori

        elif self.skipmagnetfile:  # For bilayer, need monomagdict.json
            if self.magatomlist is None:
                monomagdictpath = os.path.join(self.workdir, "monomagdict.json")
                assert os.path.exists(monomagdictpath), f"The path {monomagdictpath} does not exist."
                with open(monomagdictpath, "r") as f:
                    self.monomagdict = json.load(f)

                magnet_dict = {}
                assert os.path.exists(self.maglist), "The magnetic list file does not exist."
                with open(self.maglist, "r") as f:
                    lines = f.readlines()
                for line in lines:
                    material = line.strip()
                    magnet_dict.update(self.skipmagnetfile_getmagnet(material, self.natomlayer1))
            else:
                magnet_dict = {}
                assert os.path.exists(self.maglist), "The magnetic list file does not exist."
                with open(self.maglist, "r") as f:
                    lines = f.readlines()
                for line in lines:
                    material = line.strip()
                    poscar = os.path.join(self.workdir, material, self.pospath.replace("CONTCAR", "POSCAR"))
                    if not os.path.exists(poscar):
                        print(f"Warning: {poscar} does not exist.")
                        continue
                    with open(poscar, "r") as f:
                        poslines = f.readlines()
                    elements = list(map(str, poslines[5].split()))
                    element_ori = elements[:]
                    elemnum = list(map(int, poslines[6].split()))
                    element = []
                    for i, elem in enumerate(element_ori):
                        for j in range(elemnum[i]):
                            element.append(elem)
                    magnet_dict[material] = {}
                    magatom = []
                    for elem_ in element:
                        if elem_ in self.magatomlist:
                            magatom.append(elem_)
                        else:
                            magatom.append(None)
                    magnet_dict[material]["magatom"] = magatom
                    magnet_dict[material]["element"] = element_ori

        else:
            raise ValueError(
                "Please provide the magnetic moments file or the structure and outcar files, the latter requiring the pymatgen package.")

        with open(os.path.join(self.workdir, "magdict.json"), "w") as f:
            json.dump(magnet_dict, f)
        return magnet_dict

    def skipmagnetfile_getmagnet(self, material, natomlayer1):
        magnet_dict = {}
        poscar = os.path.join(self.workdir, material, self.pospath.replace("CONTCAR", "POSCAR"))
        if not os.path.exists(poscar):
            print(f"Warning: {poscar} does not exist.")
            pass
        with open(poscar, "r") as f:
            poslines = f.readlines()
        elements = list(map(str, poslines[5].split()))
        element_ori = elements[:]
        elemnum = list(map(int, poslines[6].split()))
        positions = np.array([list(map(float, line.split()[:3])) for line in poslines[8:8 + sum(elemnum)]])
        indices = np.argsort(positions[:, 2])
        sorted_elements = []
        for elem, num in zip(elements, elemnum):
            sorted_elements.extend([elem] * num)
        elements = [sorted_elements[i] for i in indices]
        bl1elem = elements[:natomlayer1]
        bl2elem = elements[natomlayer1:]
        magnet_dict[material] = {}
        magatom = [None] * len(elements)
        magatomacmono = []
        if "".join(bl1elem) in self.monomagdict:
            l1magatom = [mat for mat in self.monomagdict["".join(bl1elem)]["magatom"] if mat is not None]
            for el in bl1elem:
                if el in l1magatom:
                    magatomacmono.extend([el])
                else:
                    magatomacmono.extend([None])
        else:
            magatomacmono.extend([None] * natomlayer1)
        if "".join(bl2elem) in self.monomagdict:
            l2magatom = [mat for mat in self.monomagdict["".join(bl2elem)]["magatom"] if mat is not None]
            for el in bl2elem:
                if el in l2magatom:
                    magatomacmono.extend([el])
                else:
                    magatomacmono.extend([None])
        else:
            magatomacmono.extend([None] * (len(elements) - natomlayer1))

        for i, ind in enumerate(indices):
            if magatomacmono[i] is not None:
                magatom[ind] = magatomacmono[i]
        magnet_dict[material]["magatom"] = magatom
        magnet_dict[material]["element"] = element_ori
        return magnet_dict

    def genInputfile(self, INCARbasic="INCAR-basic", supercellgen="pmg", vaspkitKpoints=True, kmseshScheme=2,
                     kmeshrv=0.04):
        """
        Generate supercell and associated INCAR files.

        INCARbasic: str
            The basic INCAR file, which should contain the "INSPIN = 2" tag.
        supercellgen: str
            The method to generate the supercell.
            "pmg": Use pymatgen to generate the supercell.
            "vaspkit": Use VASPKIT to generate the supercell.
            "manual": Manually generate the supercell, and the POSCAR221 file should be stored in the work directory.
        scaling_matrix: list
            The scaling matrix of the supercell.
        vaspkitKpoints: bool
            Whether to use VASPKIT to generate KPOINTS file.
            If False, the KPOINTS file should be prepared manually.
        kmseshScheme: int
            K-Mesh Scheme.
            1: Monkhorst-Pack Scheme; 2: Gamma-Centered Scheme; 3: Irreducible K-Points with Gamma Scheme.
        kmeshrv: float
            Kmesh-Resolved Value. Accuracy level: Gamma-Only: 0; Low: 0.06~0.04; Medium: 0.04~0.03; Fine: 0.02-0.01.
        """
        assert os.path.exists(INCARbasic), "The basic INCAR file does not exist."

        scaling_matrix = [2, 2, 1]
        scalnum = "4"
        scalmat = "".join([f"{x}" for x in scaling_matrix])

        for material in self.magdict.keys():
            calpath = os.path.join(self.workdir, material)
            assert os.path.exists(calpath), f"The path {calpath} does not exist."
            nmagatom = len([x for x in self.magdict[material]["magatom"] if x is not None])

            FMdir = os.path.join(calpath, "FM")
            AFM1dir = os.path.join(calpath, "AFM1")
            AFM2dir = os.path.join(calpath, "AFM2")
            AFM3dir = os.path.join(calpath, "AFM3")

            if nmagatom == 1:
                if supercellgen.lower() == "pmg" and self.pmg:
                    from pymatgen.core.structure import Structure

                    stfile = os.path.join(calpath, "CONTCAR")
                    if not os.path.exists(stfile):
                        stfile = os.path.join(calpath, "POSCAR")
                    assert os.path.exists(stfile), f"The file {stfile} does not exist."
                    save_path = os.path.join(calpath, "POSCAR" + scalmat)
                    st = Structure.from_file(stfile)
                    st.make_supercell(scaling_matrix)
                    try:
                        st.to(save_path, "poscar")
                    except:
                        st.to("poscar", save_path)
                elif supercellgen.lower() == "vaspkit":
                    if os.path.exists(os.path.join(calpath, "POSCAR")):
                        os.system(f"echo -e \"401\n1\n2 2 1\" | (cd {calpath} && vaspkit >/dev/null 2>&1)")
                    elif os.path.exists(os.path.join(calpath, "CONTCAR")):
                        os.system(f"echo -e \"401\n2\n2 2 1\" | (cd {calpath} && vaspkit >/dev/null 2>&1)")
                    else:
                        raise FileNotFoundError(f"Neither CONTCAR nor POSCAR is found in {calpath}.")
                    os.system(f"mv {os.path.join(calpath, 'SC221.vasp')} {os.path.join(calpath, 'POSCAR' + scalmat)}")
                elif supercellgen.lower() == "manual":
                    pass
                else:
                    raise ValueError("The supercell generation method is not supported.")

                magmom = {"FM": [], "AFM1": [], "AFM2": []}
                for elem in self.magdict[material]["magatom"]:
                    if elem is None:
                        for key in magmom.keys():
                            magmom[key].extend([scalnum + '*' + str(0)])
                    else:
                        if self.magmom_setting is not None and elem in self.magmom_setting.keys():
                            elem_mag = str(self.magmom_setting[elem])
                        elif elem in maginit["MAGMOM"].keys():
                            elem_mag = str(maginit["MAGMOM"][elem])
                        else:
                            elem_mag = str(3)
                        magmom["FM"].extend([scalnum + '*' + elem_mag])
                        magmom["AFM1"].extend([elem_mag + " " + elem_mag + " -" + elem_mag + " -" + elem_mag])
                        magmom["AFM2"].extend([elem_mag + " -" + elem_mag + " -" + elem_mag + " " + elem_mag])

                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[material]["magatom"] if x is not None]

                for dir_ in [FMdir, AFM1dir, AFM2dir]:
                    if not os.path.exists(dir_):
                        os.makedirs(dir_)
                    incarfile = os.path.join(dir_, "INCAR")
                    # shutil.copy(INCARbasic, incarfile)
                    # shutil.copy(save_path, os.path.join(dir_, "POSCAR"))
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(self.magdict[material]["element"], incar=incarfile, magatomlist=magelem,
                                u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(dir_, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                    else:
                        if vaspkitKpoints:
                            os.system(f"echo -e 103 | (cd {dir_} && vaspkit >/dev/null 2>&1)")

                    if vaspkitKpoints:
                        os.system(
                            f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {dir_} && vaspkit >/dev/null 2>&1)")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        if dir_ == FMdir:
                            f.write(" ".join(magmom["FM"]))
                        elif dir_ == AFM1dir:
                            f.write(" ".join(magmom["AFM1"]))
                        elif dir_ == AFM2dir:
                            f.write(" ".join(magmom["AFM2"]))
                        f.write("\n")
                # if vaspkitKpoints:
                #     os.system(f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {FMdir} && vaspkit >/dev/null 2>&1)")
                #     for dir_ in [AFM1dir, AFM2dir]:
                #         os.system(f"ln -nfs ../FM/KPOINTS {os.path.join(dir_, 'KPOINTS')}")

            elif nmagatom == 2:
                magmom = {"FM": [], "AFM1": []}
                flag = ""
                for elem in self.magdict[material]["magatom"]:
                    if elem is None:
                        for key in magmom.keys():
                            magmom[key].extend([str(0)])
                    else:
                        if self.magmom_setting is not None and elem in self.magmom_setting.keys():
                            elem_mag = str(self.magmom_setting[elem])
                        elif elem in maginit["MAGMOM"].keys():
                            elem_mag = str(maginit["MAGMOM"][elem])
                        else:
                            elem_mag = str(3)
                        magmom["FM"].extend([elem_mag])
                        magmom["AFM1"].extend([flag + elem_mag])
                        flag = "-"

                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[material]["magatom"] if x is not None]

                if self.mkfm:
                    dirlist = [FMdir, AFM1dir]
                else:
                    dirlist = [AFM1dir]

                for dir_ in dirlist:
                    if not os.path.exists(dir_):
                        os.makedirs(dir_)
                    incarfile = os.path.join(dir_, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(self.magdict[material]["element"], incar=incarfile, magatomlist=magelem,
                                u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR')} {os.path.join(dir_, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                    else:
                        if not vaspkitKpoints:
                            print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        if dir_ == FMdir:
                            f.write(" ".join(magmom["FM"]))
                        elif dir_ == AFM1dir:
                            f.write(" ".join(magmom["AFM1"]))
                        f.write("\n")

                    if os.path.exists(os.path.join(calpath, self.kpointspath)):
                        os.system(f"cp {os.path.join(calpath, self.kpointspath)} {os.path.join(dir_, 'KPOINTS')}")
                    else:
                        if vaspkitKpoints:
                            os.system(
                                f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {dir_} && vaspkit >/dev/null 2>&1)")

            elif nmagatom == 4:  # Not tested!!!
                magmom = {"FM": [], "AFM1": [], "AFM2": [], "AFM3": []}
                flag_list = [["", "", ""], ["", "-", "-"], ["-", "-", ""], ["-", "", "-"]]
                flag_i = 0
                for elem in self.magdict[material]["magatom"]:
                    if elem is None:
                        for key in magmom.keys():
                            magmom[key].extend([str(0)])
                    else:
                        if self.magmom_setting is not None and elem in self.magmom_setting.keys():
                            elem_mag = str(self.magmom_setting[elem])
                        elif elem in maginit["MAGMOM"].keys():
                            elem_mag = str(maginit["MAGMOM"][elem])
                        else:
                            elem_mag = str(3)
                        flag_i += 1
                        magmom["FM"].extend([elem_mag])
                        magmom["AFM1"].extend([flag_list[flag_i][0] + elem_mag])
                        magmom["AFM2"].extend([flag_list[flag_i][1] + elem_mag])
                        magmom["AFM3"].extend([flag_list[flag_i][2] + elem_mag])

                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[material]["magatom"] if x is not None]

                if self.mkfm:
                    dirlist = [FMdir, AFM1dir, AFM2dir, AFM3dir]
                else:
                    dirlist = [AFM1dir, AFM2dir, AFM3dir]

                for dir_ in dirlist:
                    if not os.path.exists(dir_):
                        os.makedirs(dir_)
                    incarfile = os.path.join(dir_, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(self.magdict[material]["element"], incar=incarfile, magatomlist=magelem,
                                u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR')} {os.path.join(dir_, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                    else:
                        print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        if dir_ == FMdir:
                            f.write(" ".join(magmom["FM"]))
                        elif dir_ == AFM1dir:
                            f.write(" ".join(magmom["AFM1"]))
                        elif dir_ == AFM2dir:
                            f.write(" ".join(magmom["AFM2"]))
                        elif dir_ == AFM3dir:
                            f.write(" ".join(magmom["AFM3"]))
                        f.write("\n")
                    os.system(f"ln -nfs ../{self.kpointspath} {os.path.join(dir_, 'KPOINTS')}")

            else:
                raise ValueError("The number of magnetic atoms is not supported.")

    def catfile(self):
        return self.maglist

    def genRunsh(self, pbs=True, moduleload="vasp/5.4.4-intel2019", vasp="vasp_std_2D"):
        magcatfile = self.catfile()
        cd_com = "../" * self.multilevel
        warnings.warn("Specific calculation configurations must be reset manually.")
        with open(self.runsh, "w", newline="\n") as f:
            f.write("#!/bin/bash\n")
            if pbs:
                if subdict is not None and "PBS_set" in subdict.keys():
                    for key, value in subdict["PBS_set"].items():
                        f.write(f"#{key}{value}\n")
                    f.write("cd $PBS_O_WORKDIR\n\n")
                    if "PBS_echo" in subdict.keys():
                        for echo_ in subdict["PBS_echo"]:
                            f.write(f"{echo_}\n")
                else:
                    f.write("#PBS -N vasp\n")
                    f.write("#PBS -l nodes=1:ppn=24\n")
                    f.write("#PBS -j oe\n")
                    f.write("#PBS -q short\n")
                    f.write("#define variables\n\n")
                    f.write("cd $PBS_O_WORKDIR\n\n")
            else:
                if subdict is not None and "SLURM_set" in subdict.keys():
                    for key, value in subdict["SLURM_set"].items():
                        f.write(f"#{key}{value}\n")
                    f.write("cd $SLURM_SUBMIT_DIR\n\n")
                    if "SLURM_echo" in subdict.keys():
                        for echo_ in subdict["SLURM_echo"]:
                            f.write(f"{echo_}\n")
                else:
                    f.write("#SBATCH -J vasp\n")
                    f.write("#SBATCH -N 1\n")
                    f.write("#SBATCH -n 48\n")
                    f.write("#SBATCH -p normal\n")
                    f.write("#SBATCH -o %j.log\n")
                    f.write("#SBATCH -e %j.err\n")
                    f.write("#SBATCH -t 48:00:00\n\n")
                    f.write("cd $SLURM_SUBMIT_DIR\n\n")

            # f.write("module add vasp/5.4.4-intel2019\n\n")
            if moduleload is not None:
                module_ = moduleload.split(",")
                f.write(f"module load")
                for m in module_:
                    f.write(f" {m}")
                f.write("\n")

            f.write("MPI=mpirun\n")
            f.write(f"VASP={vasp}\n")

            f.write(f"for dir in `cat {magcatfile}`; do\n")
            f.write(f"    cd $dir\n")
            if self.skipmagnetfile:
                f.write("    $MPI $VASP > runlog\n")
                f.write("    if [ `grep -c \"or copy CONTCAR\" runlog` -ne 0 ];then\n")
                f.write("        cp CONTCAR POSCAR\n")
                f.write("        $MPI $VASP > runlog\n")
                f.write("    fi\n")
            f.write(f"    folders=(\"FM\" \"AFM1\" \"AFM2\" \"AFM3\")\n")
            f.write("    for folder in \"${folders[@]}\"; do\n")
            f.write("        if [ -d $folder ]; then\n")
            f.write("            cd $folder\n")
            f.write("            $MPI $VASP > runlog\n")
            f.write("            if [ `grep -c \"or copy CONTCAR\" runlog` -ne 0 ];then\n")
            f.write("                cp CONTCAR POSCAR\n")
            f.write("                $MPI $VASP > runlog\n")
            f.write("            fi\n")
            f.write("            cd ../\n")
            f.write("        fi\n")
            f.write("    done\n")
            f.write(f"    cd {cd_com}\n")
            f.write("done\n")

    def getEnergy(self, oszipath=None, FMm4=False, copyCONTCAR=True):
        """
        Get the energy of the magnetic states.

        oszipath: str
            The path of the OSZICAR file.
            If the saveall=True calculation is performed, the OSZICAR file will be saved in a folder named after the
            task, in which case you must specify the task name, e.g. 'opt' or 'scf'.
        FMm4: bool
            If True, the energy of the FM state will be multiplied by 4.
        """
        for material in self.magdict.keys():
            magE = {}
            magstate = "FM"
            magEstable = 0
            mag = None
            calpath = os.path.join(self.workdir, material)
            assert os.path.exists(calpath), f"The path {calpath} does not exist."

            stfile = os.path.join(calpath, self.pospath)
            if not os.path.exists(stfile):
                stfile = os.path.join(calpath, self.pospath.replace("CONTCAR", "POSCAR"))
            assert os.path.exists(stfile), f"The file {stfile} does not exist."
            with open(stfile, "r") as f:
                lines = f.readlines()
            elements = list(map(str, lines[5].split()))
            elemnum = list(map(int, lines[6].split()))
            positions = np.array([list(map(float, line.split()[:3])) for line in lines[8:8 + sum(elemnum)]])
            indices = np.argsort(positions[:, 2])
            sorted_elements = []
            for elem, num in zip(elements, elemnum):
                sorted_elements.extend([elem] * num)
            elements = [sorted_elements[i] for i in indices]
            elements_inv = elements[::-1]
            mlabel = "".join(elements)
            mlabel_inv = "".join(elements_inv)
            self.magdict[material]["mlabel"] = mlabel
            self.magdict[material]["mlabel_inv"] = mlabel_inv

            FMdir = os.path.join(calpath, "FM")
            AFM1dir = os.path.join(calpath, "AFM1")
            AFM2dir = os.path.join(calpath, "AFM2")
            AFM3dir = os.path.join(calpath, "AFM3")

            incarpath = None
            if os.path.exists(FMdir) and os.path.exists(os.path.join(FMdir, "OSZICAR")):
                oszicar = os.path.join(FMdir, "OSZICAR")
                contcar = os.path.join(FMdir, "CONTCAR")
                if os.path.exists(oszicar):
                    grepE = subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $5}'",
                                           shell=True, capture_output=True, text=True)
                    magE["FM"] = float(grepE.stdout.strip())
                    if magE["FM"] < magEstable:
                        magEstable = magE["FM"]
                        mag = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $10}'",
                                                   shell=True, capture_output=True, text=True).stdout.strip())
            else:
                if oszipath:
                    oszicar = os.path.join(calpath, oszipath, "OSZICAR")
                    contcar = os.path.join(calpath, oszipath, "CONTCAR")
                else:
                    oszicar = os.path.join(calpath, "OSZICAR")
                    contcar = os.path.join(calpath, "CONTCAR")
                if os.path.exists(oszicar):
                    grepE = subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $5}'",
                                           shell=True, capture_output=True, text=True)
                    magE["FM"] = float(grepE.stdout.strip()) if FMm4 is False else float(grepE.stdout.strip()) * 4
                    if magE["FM"] < magEstable:
                        magEstable = magE["FM"]
                        mag = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $10}'",
                                                   shell=True, capture_output=True, text=True).stdout.strip())

            outcarpath = os.path.join(os.path.dirname(oszicar), "OUTCAR")
            dirlist = ["AFM1", "AFM2", "AFM3"]
            for i, dir_ in enumerate([AFM1dir, AFM2dir, AFM3dir]):
                if os.path.exists(dir_):
                    oszicar = os.path.join(dir_, "OSZICAR")
                    if os.path.exists(oszicar):
                        grepE = subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $5}'",
                                               shell=True, capture_output=True, text=True)
                        magE[dirlist[i]] = float(grepE.stdout.strip())
                        if magE[dirlist[i]] < magEstable:
                            magEstable = magE[dirlist[i]]
                            magstate = dirlist[i]
                            contcar = os.path.join(dir_, "CONTCAR")
                            incarpath = os.path.join(dir_, "INCAR")
                            outcarpath = os.path.join(dir_, "OUTCAR")
                            mag = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $10}'",
                                                       shell=True, capture_output=True, text=True).stdout.strip())

            self.magdict[material]["magE"] = magE
            self.magdict[material]["magstate"] = magstate
            self.magdict[material]["mag"] = mag

            if outcarpath and os.path.exists(outcarpath):
                if self.Outcar:
                    magmom = self.Outcar(outcarpath).magnetization
                    magmom_list = [magtot['tot'] for magtot in magmom]
                else:
                    startnu = int(subprocess.run(f"grep -n '# of ion' {outcarpath} | tail -n 1 | cut -d':' -f1",
                                                 shell=True, capture_output=True, text=True).stdout.strip())
                    endnu = int(subprocess.run(f"grep -n 'tot ' {outcarpath} | tail -n 1 | cut -d':' -f1",
                                               shell=True, capture_output=True, text=True).stdout.strip())
                    nbegin = startnu + 2
                    nend = endnu - 1
                    magmom_list = []
                    for i in range(nbegin, nend):
                        grepmagmom = subprocess.run(f"sed -n '{i}p' {outcarpath}",
                                                    shell=True, capture_output=True, text=True).stdout.strip()
                        magmom_list.append(float(grepmagmom.split()[-1]))
                magcor = determine_mag_state(magmom_list, self.maglowerl)
                self.magdict[material]["magstatecorrect"] = magcor
                if magcor == "FM" or magcor == "NM":
                    if os.path.exists(os.path.join(calpath, 'MAGMOM')):
                        os.system(f"rm {os.path.join(calpath, 'MAGMOM')}")
                    incarpath = None

            if incarpath:
                grepMAG = subprocess.run(f"grep 'MAGMOM' {incarpath}",
                                         shell=True, capture_output=True, text=True)
                magcom = grepMAG.stdout.strip()
                if magcom:
                    with open(os.path.join(calpath, 'MAGMOM'), "w", newline="\n") as f:
                        f.write(f"{magcom}")

            if copyCONTCAR and contcar:
                os.system(f"cp {contcar} {os.path.join(calpath, 'POSCAR')}")

        magdictsave = {}
        for material in self.magdict.keys():
            mlabel_ = self.magdict[material]["mlabel"]
            mlabel_inv_ = self.magdict[material]["mlabel_inv"]
            magdictsave[mlabel_] = {}
            magdictsave[mlabel_]["material"] = material
            magdictsave[mlabel_]["magatom"] = self.magdict[material]["magatom"]
            magdictsave[mlabel_]["magE"] = self.magdict[material]["magE"]
            magdictsave[mlabel_]["magstate"] = self.magdict[material]["magstate"]
            magdictsave[mlabel_]["mag"] = self.magdict[material]["mag"]
            if "magstatecorrect" in self.magdict[material]:
                magdictsave[mlabel_]["magstatecorrect"] = self.magdict[material]["magstatecorrect"]
            if mlabel_inv_ != mlabel_:
                magdictsave[mlabel_inv_] = magdictsave[mlabel_]
        with open(os.path.join(self.workdir, "monomagdict.json"), "w") as f:
            json.dump(magdictsave, f)


class BiMagState(MonoMagState):
    """
    The number of layers is  2. Magnetic ground state determination is performed in layers, with the monolayer
    magnetic ground state first confirmed before the interlayer magnetic coupling state is determined. In this case,
    the monolayer's magnetic ground state must be provided beforehand.
    """

    def genbidict(self, mlist, natomlayer1=3):
        biInputdict = {}
        biwarn = []
        biwarnmid = []
        for material in mlist:
            magatom = self.magdict[material]["magatom"]
            materialsplit = material.split("/")[0].split("_")
            mid = materialsplit[0]
            self.magdict[material]["mid"] = mid
            if mid not in biInputdict.keys():
                biInputdict[mid] = {}
                biInputdict[mid]["Inputpath"] = material
                biInputdict[mid]["magatom"] = magatom
                biInputdict[mid]["submaterial"] = [material]
                natomlayer1 = int(materialsplit[2]) if len(materialsplit) == 3 else natomlayer1
                if natomlayer1 != self.natomlayer1:
                    magatom = self.skipmagnetfile_getmagnet(material, natomlayer1)[material]["magatom"]
                    biInputdict[mid]["magatom"] = magatom
                magmomlist = []
                for i, elem in enumerate(magatom):
                    if elem is not None:
                        if self.magmom_setting is not None and elem in self.magmom_setting.keys():
                            magmomlist.append(self.magmom_setting[elem])
                        elif elem in maginit["MAGMOM"].keys():
                            magmomlist.append(maginit["MAGMOM"][elem])
                        else:
                            magmomlist.append(3)
                    else:
                        magmomlist.append(0)
                biInputdict[mid]["natomlayer1"] = natomlayer1
                biInputdict[mid]["magmomlist"] = magmomlist

                calpath = os.path.join(self.workdir, material)
                stfile = os.path.join(calpath, self.pospath)
                if not os.path.exists(stfile):
                    stf_ = self.pospath.replace("CONTCAR", "POSCAR")
                    stfile = os.path.join(calpath, stf_)
                assert os.path.exists(stfile), f"The file {stfile} does not exist."
                with open(stfile, "r") as f:
                    lines = f.readlines()
                elements = list(map(str, lines[5].split()))
                element_ori = elements[:]
                elemnum = list(map(int, lines[6].split()))
                positions = np.array([list(map(float, line.split()[:3])) for line in lines[8:8 + sum(elemnum)]])
                indices = np.argsort(positions[:, 2])
                sorted_elements = []
                for elem, num in zip(elements, elemnum):
                    sorted_elements.extend([elem] * num)
                elements = [sorted_elements[i] for i in indices]
                magmomlistsorted = [magmomlist[i] for i in indices]
                blabel = "".join(elements)
                bl1label = "".join(elements[:natomlayer1])
                bl2label = "".join(elements[natomlayer1:])
                biInputdict[mid]["blabel"] = blabel
                biInputdict[mid]["bl1label"] = bl1label
                biInputdict[mid]["bl2label"] = bl2label
                biInputdict[mid]["element"] = element_ori

                layer1nmag = len([x for x in magmomlistsorted[:natomlayer1] if x != 0])
                layer2nmag = len([x for x in magmomlistsorted[natomlayer1:] if x != 0])
                biInputdict[mid]["nmag"] = layer1nmag + layer2nmag

                def getmagstate(nmag, label):
                    if nmag == 0:
                        return 0  # non-magnetic
                    else:
                        if label in self.monomagdict.keys():
                            magstate_ = self.monomagdict[label]["magstate"]
                            if magstate_ == "FM":
                                return 3
                            elif magstate_ == "AFM1":
                                return 1
                            elif magstate_ == "AFM2":
                                return 2
                            else:
                                return 4  # AFM3
                        else:
                            # warnings.warn(f"The monolayer magnetic ground state of {label} is not provided.")
                            return 5  # unknown

                magstate1 = getmagstate(layer1nmag, bl1label)
                magstate2 = getmagstate(layer2nmag, bl2label)
                magstate = str(magstate1) + str(magstate2) if magstate1 < magstate2 else str(magstate2) + str(magstate1)
                biInputdict[mid]["magstate"] = magstate
                ranking = np.empty_like(indices)
                ranking[indices] = np.arange(len(indices))
                biInputdict[mid]["atommagstate"] = [magstate1 if i_ < natomlayer1 else magstate2 for i_ in ranking]
                biInputdict[mid]["atomlayer"] = [1 if i_ < natomlayer1 else 2 for i_ in ranking]
            else:
                if not self.skipmagnetfile:
                    magatom_bi = biInputdict[mid]["magatom"]
                    if magatom != magatom_bi:
                        # warnings.warn(f"The magnetic atom {magatom} in {material} is different from the magnetic atom "
                        #               f"{magatom_bi} in {biInputdict[mid]['Inputpath']}.")
                        biwarn.append(material)
                        biwarnmid.append(mid)
                    else:
                        biInputdict[mid]["submaterial"].append(material)
                else:
                    biInputdict[mid]["submaterial"].append(material)
        return biInputdict, biwarn, biwarnmid

    def iniMonoMagstates(self, natomlayer1=3):
        self.bimaglist = "bimaglist.txt"

        if not os.path.exists(os.path.join(self.workdir, 'biInputdict.json')):
            if not self.skipmagnetfile and not os.path.exists(os.path.join(self.workdir, 'monomagdict.json')):
                monomagdictpath = os.path.join(self.workdir, "monomagdict.json")
                assert os.path.exists(monomagdictpath), f"The path {monomagdictpath} does not exist."
                with open(monomagdictpath, "r") as f:
                    self.monomagdict = json.load(f)
            self.biInputdict, biwarn, biwarnmid = self.genbidict(self.magdict.keys(), natomlayer1)
            i_ = 0
            while len(biwarn) > 0:
                i_ += 1
                self.biInputdict2, biwarn2, biwarnmid2 = self.genbidict(biwarn, natomlayer1)
                for mid in self.biInputdict2.keys():
                    if mid in self.biInputdict.keys():
                        self.biInputdict[mid + "_" + str(i_)] = self.biInputdict2[mid]
                biwarn = biwarn2[:]
            with open(os.path.join(self.workdir, 'biInputdict.json'), 'w') as f:
                json.dump(self.biInputdict, f)
        else:
            with open(os.path.join(self.workdir, 'biInputdict.json'), 'r') as f:
                self.biInputdict = json.load(f)

    def gensupercell(self, pathlist, scaling_matrix, scalmat, supercellgen="pmg"):
        if supercellgen.lower() == "pmg":
            from pymatgen.core.structure import Structure

            for path_ in pathlist:
                calpath = os.path.join(self.workdir, path_)
                stfile = os.path.join(calpath, "CONTCAR")
                if not os.path.exists(stfile):
                    stfile = os.path.join(calpath, "POSCAR")
                assert os.path.exists(stfile), f"The file {stfile} does not exist."
                save_path = os.path.join(calpath, "POSCAR" + scalmat)
                st = Structure.from_file(stfile)
                st.make_supercell(scaling_matrix)
                try:
                    st.to(save_path, "poscar")
                except:
                    st.to("poscar", save_path)
        elif supercellgen.lower() == "vaspkit":
            for path_ in pathlist:
                calpath = os.path.join(self.workdir, path_)
                if os.path.exists(os.path.join(calpath, "POSCAR")):
                    os.system(f"echo -e \"401\n1\n2 2 1\" | (cd {calpath} && vaspkit >/dev/null 2>&1)")
                elif os.path.exists(os.path.join(calpath, "CONTCAR")):
                    os.system(f"echo -e \"401\n2\n2 2 1\" | (cd {calpath} && vaspkit >/dev/null 2>&1)")
                else:
                    raise FileNotFoundError(f"Neither CONTCAR nor POSCAR is found in {calpath}.")
                os.system(
                    f"mv {os.path.join(calpath, 'SC221.vasp')} {os.path.join(calpath, 'POSCAR' + scalmat)}")
        elif supercellgen.lower() == "manual":
            pass
        else:
            raise ValueError("The supercell generation method is not supported.")

    def subgenInputfile(self, pathlist, dirprefix, scalmat):
        dirdict = {}

        for dir_ in dirprefix:
            dirdict[dir_] = []
            for path_ in pathlist:
                dirdict[dir_].append(os.path.join(self.workdir, path_, dir_))

        incarfile = os.path.join(dirdict[dirprefix[0]][0], "INCAR")
        kpointsfile = os.path.join(dirdict[dirprefix[0]][0], "KPOINTS")
        potcarfile = os.path.join(dirdict[dirprefix[0]][0], "POTCAR")

        for file_ in [incarfile, kpointsfile, potcarfile]:
            if not os.path.exists(file_):
                warnings.warn(f"The file {file_} does not exist.")

        for i_, path_ in enumerate(dirdict[dirprefix[0]][1:]):
            if not os.path.exists(path_):
                os.makedirs(path_)
            incarsave = os.path.join(path_, "INCAR")
            os.system(f"cp {incarfile} {incarsave}")
            kpointssave = os.path.join(path_, "KPOINTS")
            os.system(f"cp {kpointsfile} {kpointssave}")
            potcarsave = os.path.join(path_, "POTCAR")
            os.system(f"cp {potcarfile} {potcarsave}")
            os.system(f"cp {os.path.join(pathlist[i_ + 1], 'POSCAR' + scalmat)} {os.path.join(path_, 'POSCAR')}")

        for dir_ in dirprefix[1:]:
            incarfile = os.path.join(dirdict[dir_][0], "INCAR")
            for i_, path_ in enumerate(dirdict[dir_][1:]):
                if not os.path.exists(path_):
                    os.makedirs(path_)
                incarsave = os.path.join(path_, "INCAR")
                os.system(f"cp {incarfile} {incarsave}")
                os.system(f"ln -nfs ../{dirprefix[0]}/KPOINTS {os.path.join(path_, 'KPOINTS')}")
                os.system(f"ln -nfs ../{dirprefix[0]}/POTCAR {os.path.join(path_, 'POTCAR')}")
                os.system(f"cp {os.path.join(pathlist[i_ + 1], 'POSCAR' + scalmat)} {os.path.join(path_, 'POSCAR')}")

    def genInputfile(self, INCARbasic="INCAR-basic", supercellgen="pmg", vaspkitKpoints=True, kmseshScheme=2,
                     kmeshrv=0.04):
        assert os.path.exists(INCARbasic), "The basic INCAR file does not exist."

        scaling_matrix = [2, 2, 1]
        scalnum = "4"
        scalmat = "".join([f"{x}" for x in scaling_matrix])

        for mid in self.biInputdict.keys():
            magstate = self.biInputdict[mid]["magstate"]
            nmag = self.biInputdict[mid]["nmag"]
            calpath = os.path.join(self.workdir, self.biInputdict[mid]["Inputpath"])
            pathlist = self.biInputdict[mid]["submaterial"]
            element = self.biInputdict[mid]["element"]

            if magstate == "03" or (magstate == "00" and nmag == 0):
                continue

            elif magstate == "01" or magstate == "02":
                AFMdir = {"AFM1": os.path.join(calpath, "AFM1"), "AFM2": os.path.join(calpath, "AFM2")}
                if nmag == 1:
                    self.bicalmaterial.extend(pathlist)
                    if magstate == "01":
                        dirindex = "AFM1"
                    else:
                        dirindex = "AFM2"
                    self.gensupercell(pathlist, scaling_matrix, scalmat, supercellgen)
                    magmom = []
                    for elemmag in self.biInputdict[mid]["magmomlist"]:
                        if elemmag == 0:
                            magmom.extend([scalnum + '*' + str(0)])
                        else:
                            if dirindex == "AFM1":
                                magmom.extend(
                                    [str(elemmag) + " " + str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag)])
                            else:
                                magmom.extend(
                                    [str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag) + " " + str(elemmag)])
                    magelem = None
                    if self.onlymagpu:
                        magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                    if not os.path.exists(AFMdir[dirindex]):
                        os.makedirs(AFMdir[dirindex])
                    incarfile = os.path.join(AFMdir[dirindex], "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(
                        f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(AFMdir[dirindex], 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(AFMdir[dirindex], 'POTCAR')}")
                    else:
                        if not vaspkitKpoints:
                            print(f"The POTCAR file of {calpath} does not exist.")
                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        f.write(" ".join(magmom))
                        f.write("\n")
                    if vaspkitKpoints:
                        os.system(
                            f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {AFMdir[dirindex]} && vaspkit >/dev/null 2>&1)")
                    self.subgenInputfile(pathlist, [dirindex], scalmat)
                elif nmag == 2:
                    self.bicalmaterial.extend(pathlist)
                    magmom = []
                    flag = ""
                    for elemmag in self.biInputdict[mid]["magmomlist"]:
                        if elemmag == 0:
                            magmom.extend([str(0)])
                        else:
                            elemmag = flag + str(elemmag)
                            flag = "-"
                            magmom.extend([elemmag])
                    magelem = None
                    if self.onlymagpu:
                        magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                    if not os.path.exists(AFMdir['AFM1']):
                        os.makedirs(AFMdir['AFM1'])
                    incarfile = os.path.join(AFMdir['AFM1'], "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR')} {os.path.join(AFMdir['AFM1'], 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(AFMdir['AFM1'], 'POTCAR')}")
                    else:
                        if not vaspkitKpoints:
                            print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        f.write(" ".join(magmom))
                        f.write("\n")
                    os.system(f"ln -nfs ../{self.kpointspath} {os.path.join(AFMdir['AFM1'], 'KPOINTS')}")
                    self.subgenInputfile(pathlist, ["AFM1"], "")
                else:
                    print(f"Systems with more than 2 magnetic atoms are not currently supported.")
                    print(f"Please check the magnetic atom list of {mid}.")

            elif magstate == "05":
                FMdir = os.path.join(calpath, "FM")
                AFM1dir = os.path.join(calpath, "AFM1")
                AFM2dir = os.path.join(calpath, "AFM2")
                if nmag == 1:
                    self.bicalmaterial.extend(pathlist)
                    self.gensupercell(pathlist, scaling_matrix, scalmat, supercellgen)
                    magmom = {"FM": [], "AFM1": [], "AFM2": []}
                    for elemmag in self.biInputdict[mid]["magmomlist"]:
                        if elemmag == 0:
                            for key in magmom.keys():
                                magmom[key].extend([scalnum + '*' + str(0)])
                        else:
                            magmom["FM"].extend([scalnum + '*' + str(elemmag)])
                            magmom["AFM1"].extend(
                                [str(elemmag) + " " + str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag)])
                            magmom["AFM2"].extend(
                                [str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag) + " " + str(elemmag)])
                    magelem = None
                    if self.onlymagpu:
                        magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                    for dir_ in [FMdir, AFM1dir, AFM2dir]:
                        if not os.path.exists(dir_):
                            os.makedirs(dir_)
                        incarfile = os.path.join(dir_, "INCAR")
                        os.system(f"cp {INCARbasic} {incarfile}")
                        if self.pU:
                            INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                    openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                        os.system(f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(dir_, 'POSCAR')}")
                        if os.path.exists(os.path.join(calpath, "POTCAR")):
                            os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                        else:
                            if not vaspkitKpoints:
                                print(f"The POTCAR file of {calpath} does not exist.")

                        with open(incarfile, "a") as f:
                            f.write("MAGMOM = ")
                            if dir_ == FMdir:
                                f.write(" ".join(magmom["FM"]))
                            elif dir_ == AFM1dir:
                                f.write(" ".join(magmom["AFM1"]))
                            elif dir_ == AFM2dir:
                                f.write(" ".join(magmom["AFM2"]))
                            f.write("\n")
                    if vaspkitKpoints:
                        os.system(
                            f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {FMdir} && vaspkit >/dev/null 2>&1)")
                        for dir_ in [AFM1dir, AFM2dir]:
                            os.system(f"ln -nfs ../FM/KPOINTS {os.path.join(dir_, 'KPOINTS')}")
                    self.subgenInputfile(pathlist, ["FM", "AFM1", "AFM2"], scalmat)
                elif nmag == 2:
                    self.bicalmaterial.extend(pathlist)
                    magmom = []
                    flag = ""
                    for elemmag in self.biInputdict[mid]["magmomlist"]:
                        if elemmag == 0:
                            magmom.extend([str(0)])
                        else:
                            elemmag = flag + str(elemmag)
                            flag = "-"
                            magmom.extend([elemmag])

                    magelem = None
                    if self.onlymagpu:
                        magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                    if not os.path.exists(AFM1dir):
                        os.makedirs(AFM1dir)
                    incarfile = os.path.join(AFM1dir, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR')} {os.path.join(AFM1dir, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(AFM1dir, 'POTCAR')}")
                    else:
                        print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        f.write(" ".join(magmom))
                        f.write("\n")
                    os.system(f"ln -nfs ../{self.kpointspath} {os.path.join(AFM1dir, 'KPOINTS')}")
                    self.subgenInputfile(pathlist, ["AFM1"], "")
                else:
                    print(f"Systems with more than 2 magnetic atoms are not currently supported.")
                    print(f"Please check the magnetic atom list of {mid}.")

            elif magstate == "35":
                AFM1dir = os.path.join(calpath, "AFM1")
                magmomlist = self.biInputdict[mid]["magmomlist"]
                atomlayer = self.biInputdict[mid]["atomlayer"]
                atommagstate = self.biInputdict[mid]["atommagstate"]
                l1magstate = [x for i, x in enumerate(atommagstate) if atomlayer[i] == 1][0]
                if l1magstate == 5:
                    unkownnmag = len([x for i, x in enumerate(magmomlist) if atomlayer[i] == 1 and x != 0])
                else:
                    unkownnmag = len([x for i, x in enumerate(magmomlist) if atomlayer[i] == 2 and x != 0])
                if unkownnmag <= 2:
                    self.bicalmaterial.extend(pathlist)
                    magmom = []
                    flag = ""
                    for i_, elemmag in enumerate(magmomlist):
                        if elemmag == 0:
                            magmom.extend([str(0)])
                        else:
                            if unkownnmag == 1:
                                if atommagstate[i_] == 5:
                                    elemmag = "-" + str(elemmag)
                                else:
                                    elemmag = str(elemmag)
                                magmom.extend([elemmag])
                            else:
                                if atommagstate[i_] == 5:
                                    elemmag = flag + str(elemmag)
                                    flag = "-"
                                else:
                                    elemmag = str(elemmag)
                                magmom.extend([elemmag])
                    magelem = None
                    if self.onlymagpu:
                        magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                    if not os.path.exists(AFM1dir):
                        os.makedirs(AFM1dir)
                    incarfile = os.path.join(AFM1dir, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR')} {os.path.join(AFM1dir, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(AFM1dir, 'POTCAR')}")
                    else:
                        print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        f.write(" ".join(magmom))
                        f.write("\n")
                    os.system(f"ln -nfs ../{self.kpointspath} {os.path.join(AFM1dir, 'KPOINTS')}")
                    self.subgenInputfile(pathlist, ["AFM1"], "")
                else:
                    print(f"Systems with more than 2 magnetic atoms are not currently supported.")
                    print(f"Please check the magnetic atom list of {mid}.")

            elif magstate == "55":
                AFM1dir = os.path.join(calpath, "AFM1")
                magmomlist = self.biInputdict[mid]["magmomlist"]
                atomlayer = self.biInputdict[mid]["atomlayer"]
                l1nmag = len([x for i, x in enumerate(magmomlist) if atomlayer[i] == 1 and x != 0])
                if l1nmag == 2:
                    layerindex = 1
                else:
                    layerindex = 2
                if nmag == 2 or nmag == 3:
                    self.bicalmaterial.extend(pathlist)
                    magmom = []
                    flag = ""
                    for i_, elemmag in enumerate(magmomlist):
                        if elemmag == 0:
                            magmom.extend([str(0)])
                        else:
                            if nmag == 2:
                                if atomlayer[i_] == layerindex:
                                    elemmag = "-" + str(elemmag)
                                else:
                                    elemmag = str(elemmag)
                                magmom.extend([elemmag])
                            else:
                                if atomlayer[i_] == layerindex:
                                    elemmag = flag + str(elemmag)
                                    flag = "-"
                                else:
                                    elemmag = str(elemmag)
                                magmom.extend([elemmag])
                    magelem = None
                    if self.onlymagpu:
                        magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                    if not os.path.exists(AFM1dir):
                        os.makedirs(AFM1dir)
                    incarfile = os.path.join(AFM1dir, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR')} {os.path.join(AFM1dir, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(AFM1dir, 'POTCAR')}")
                    else:
                        print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        f.write(" ".join(magmom))
                        f.write("\n")
                    os.system(f"ln -nfs ../{self.kpointspath} {os.path.join(AFM1dir, 'KPOINTS')}")
                    self.subgenInputfile(pathlist, ["AFM1"], "")
                elif nmag == 4:
                    AFM2dir = os.path.join(calpath, "AFM2")
                    AFM3dir = os.path.join(calpath, "AFM3")
                    magmom = {"AFM1": [], "AFM2": [], "AFM3": []}

                    flag_list_l1 = ["", "", ""]
                    flag_list_l2 = ["", "", ""]
                    for i_, elemmag in enumerate(magmomlist):
                        if elemmag == 0:
                            for key in magmom.keys():
                                magmom[key].extend([str(0)])
                        else:
                            if atomlayer[i_] == 1:
                                flag_list = flag_list_l1[:]
                                flag_list_l1 = ["", "-", "-"]
                            else:
                                flag_list = flag_list_l2[:]
                                flag_list_l2 = ["-", "", "-"]
                            magmom["AFM3"].extend([flag_list[0] + str(elemmag)])
                            magmom["AFM2"].extend([flag_list[1] + str(elemmag)])
                            magmom["AFM1"].extend([flag_list[2] + str(elemmag)])
                    magelem = None
                    if self.onlymagpu:
                        magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                    for dir_ in [AFM1dir, AFM2dir, AFM3dir]:
                        if not os.path.exists(dir_):
                            os.makedirs(dir_)
                        incarfile = os.path.join(dir_, "INCAR")
                        os.system(f"cp {INCARbasic} {incarfile}")
                        if self.pU:
                            INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                    openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                        os.system(f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(dir_, 'POSCAR')}")
                        if os.path.exists(os.path.join(calpath, "POTCAR")):
                            os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                        else:
                            print(f"The POTCAR file of {calpath} does not exist.")
                        with open(incarfile, "a") as f:
                            f.write("MAGMOM = ")
                            if dir_ == AFM1dir:
                                f.write(" ".join(magmom["AFM1"]))
                            elif dir_ == AFM2dir:
                                f.write(" ".join(magmom["AFM2"]))
                            elif dir_ == AFM3dir:
                                f.write(" ".join(magmom["AFM3"]))
                            f.write("\n")
                        os.system(f"ln -nfs ../{self.kpointspath} {os.path.join(dir_, 'KPOINTS')}")
                    self.subgenInputfile(pathlist, ["AFM1", "AFM2", "AFM3"], "")
                else:
                    print(f"Systems with more than 2 magnetic atoms are not currently supported.")
                    print(f"Please check the magnetic atom list of {mid}.")

            elif magstate == "33":
                self.bicalmaterial.extend(pathlist)
                FMdir = None
                if not os.path.exists(os.path.join(calpath, "FM")) and not os.path.exists(
                        os.path.join(calpath, self.outpath)):
                    FMdir = os.path.join(calpath, "FM")
                AFM1dir = os.path.join(calpath, "AFM1")
                magmom = {"FM": [], "AFM1": []}
                for i_, elemmag in enumerate(self.biInputdict[mid]["magmomlist"]):
                    if elemmag == 0:
                        for key in magmom.keys():
                            magmom[key].extend([str(0)])
                    else:
                        if self.biInputdict[mid]["atomlayer"][i_] == 1:
                            flag = ""
                        else:
                            flag = "-"
                        magmom["FM"].extend([str(elemmag)])
                        magmom["AFM1"].extend([flag + str(elemmag)])

                dirlist = []
                if FMdir or self.mkfm:
                    if not os.path.exists(FMdir):
                        os.makedirs(FMdir)
                    dirlist.append(FMdir)
                if not os.path.exists(AFM1dir):
                    os.makedirs(AFM1dir)
                dirlist.append(AFM1dir)
                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                for dir_ in dirlist:
                    incarfile = os.path.join(dir_, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR')} {os.path.join(dir_, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                    elif vaspkitKpoints:
                        os.system(f'echo -e "103" | (cd {dir_} && vaspkit >/dev/null 2>&1)')
                    else:
                        print(f"The POTCAR file of {calpath} does not exist.")
                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        f.write(" ".join(magmom["FM"]))
                        f.write("\n")
                    if os.path.exists(os.path.join(os.path.dirname(dir_), self.kpointspath)):
                        os.system(f"ln -nfs ../KPOINTS {os.path.join(dir_, self.kpointspath)}")
                    elif vaspkitKpoints:
                        os.system(
                            f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {dir_} && vaspkit >/dev/null 2>&1)")
                    else:
                        print(f"The KPOINTS file of {calpath} does not exist.")
                if FMdir:
                    self.subgenInputfile(pathlist, ["FM", "AFM1"], "")
                else:
                    self.subgenInputfile(pathlist, ["AFM1"], "")

            elif (magstate == "13" or magstate == "23") and nmag == 2:
                self.bicalmaterial.extend(pathlist)
                FMdir = os.path.join(calpath, "FM")
                self.gensupercell(pathlist, scaling_matrix, scalmat, supercellgen)
                magmom = []
                for i_, elemmag in enumerate(self.biInputdict[mid]["magmomlist"]):
                    if elemmag == 0:
                        magmom.extend([scalnum + '*' + str(0)])
                    else:
                        if self.biInputdict[mid]["atommagstate"][i_] == 1:
                            magmom.extend(
                                [str(elemmag) + " " + str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag)])
                        elif self.biInputdict[mid]["atommagstate"][i_] == 2:
                            magmom.extend(
                                [str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag) + " " + str(elemmag)])
                        else:
                            magmom.extend([scalnum + '*' + str(elemmag)])
                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                if not os.path.exists(FMdir):
                    os.makedirs(FMdir)
                incarfile = os.path.join(FMdir, "INCAR")
                os.system(f"cp {INCARbasic} {incarfile}")
                if self.pU:
                    INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                            openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                os.system(f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(FMdir, 'POSCAR')}")
                if os.path.exists(os.path.join(calpath, "POTCAR")):
                    os.system(f"ln -nfs ../POTCAR {os.path.join(FMdir, 'POTCAR')}")
                else:
                    print(f"The POTCAR file of {calpath} does not exist.")

                with open(incarfile, "a") as f:
                    f.write("MAGMOM = ")
                    f.write(" ".join(magmom))
                    f.write("\n")
                if vaspkitKpoints:
                    os.system(
                        f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {FMdir} && vaspkit >/dev/null 2>&1)")
                self.subgenInputfile(pathlist, ["FM"], scalmat)

            elif magstate == "12" and nmag == 2:
                self.bicalmaterial.extend(pathlist)
                AFM1dir = os.path.join(calpath, "AFM1")
                self.gensupercell(pathlist, scaling_matrix, scalmat, supercellgen)
                magmom = []
                for i_, elemmag in enumerate(self.biInputdict[mid]["magmomlist"]):
                    if elemmag == 0:
                        magmom.extend([scalnum + '*' + str(0)])
                    else:
                        if self.biInputdict[mid]["atommagstate"][i_] == 1:
                            magmom.extend(
                                [str(elemmag) + " " + str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag)])
                        else:
                            magmom.extend(
                                [str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag) + " " + str(elemmag)])
                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                if not os.path.exists(AFM1dir):
                    os.makedirs(AFM1dir)
                incarfile = os.path.join(AFM1dir, "INCAR")
                os.system(f"cp {INCARbasic} {incarfile}")
                if self.pU:
                    INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                            openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                os.system(f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(AFM1dir, 'POSCAR')}")
                if os.path.exists(os.path.join(calpath, "POTCAR")):
                    os.system(f"ln -nfs ../POTCAR {os.path.join(AFM1dir, 'POTCAR')}")
                else:
                    print(f"The POTCAR file of {calpath} does not exist.")

                with open(incarfile, "a") as f:
                    f.write("MAGMOM = ")
                    f.write(" ".join(magmom))
                    f.write("\n")
                if vaspkitKpoints:
                    os.system(
                        f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {AFM1dir} && vaspkit >/dev/null 2>&1)")
                self.subgenInputfile(pathlist, ["AFM1"], scalmat)

            elif magstate == "11" and nmag == 2:
                self.bicalmaterial.extend(pathlist)
                FMdir = os.path.join(calpath, "FM")
                AFM1dir = os.path.join(calpath, "AFM1")
                AFM2dir = os.path.join(calpath, "AFM2")
                self.gensupercell(pathlist, scaling_matrix, scalmat, supercellgen)
                magmom = {"FM": [], "AFM1": [], "AFM2": []}
                magatom_i = 0
                for elemmag in self.biInputdict[mid]["magmomlist"]:
                    if elemmag == 0:
                        for key in magmom.keys():
                            magmom[key].extend([scalnum + '*' + str(0)])
                    else:
                        if magatom_i == 0:
                            for key in magmom.keys():
                                magmom[key].extend(
                                    [str(elemmag) + " -" + str(elemmag) + " " + str(elemmag) + " -" + str(elemmag)])
                        else:
                            magmom["FM"].extend(
                                [str(elemmag) + " -" + str(elemmag) + " " + str(elemmag) + " -" + str(elemmag)])
                            magmom["AFM1"].extend(
                                ["-" + str(elemmag) + " " + str(elemmag) + " -" + str(elemmag) + " " + str(elemmag)])
                            magmom["AFM2"].extend(
                                [str(elemmag) + " " + str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag)])
                        magatom_i += 1
                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                for dir_ in [FMdir, AFM1dir, AFM2dir]:
                    if not os.path.exists(dir_):
                        os.makedirs(dir_)
                    incarfile = os.path.join(dir_, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(dir_, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                    else:
                        print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        if dir_ == FMdir:
                            f.write(" ".join(magmom["FM"]))
                        elif dir_ == AFM1dir:
                            f.write(" ".join(magmom["AFM1"]))
                        elif dir_ == AFM2dir:
                            f.write(" ".join(magmom["AFM2"]))
                        f.write("\n")
                if vaspkitKpoints:
                    os.system(
                        f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {FMdir} && vaspkit >/dev/null 2>&1)")
                    for dir_ in [AFM1dir, AFM2dir]:
                        os.system(f"ln -nfs ../FM/KPOINTS {os.path.join(dir_, 'KPOINTS')}")
                self.subgenInputfile(pathlist, ["FM", "AFM1", "AFM2"], scalmat)

            elif magstate == "22" and nmag == 2:
                self.bicalmaterial.extend(pathlist)
                FMdir = os.path.join(calpath, "FM")
                AFM1dir = os.path.join(calpath, "AFM1")
                self.gensupercell(pathlist, scaling_matrix, scalmat, supercellgen)
                magmom = {"FM": [], "AFM1": []}
                magatom_i = 0
                for elemmag in self.biInputdict[mid]["magmomlist"]:
                    if elemmag == 0:
                        for key in magmom.keys():
                            magmom[key].extend([scalnum + '*' + str(0)])
                    else:
                        if magatom_i == 0:
                            for key in magmom.keys():
                                magmom[key].extend(
                                    [str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag) + " " + str(elemmag)])
                        else:
                            magmom["FM"].extend(
                                [str(elemmag) + " -" + str(elemmag) + " -" + str(elemmag) + " " + str(elemmag)])
                            magmom["AFM1"].extend(
                                ["-" + str(elemmag) + " " + str(elemmag) + " " + str(elemmag) + " -" + str(elemmag)])
                        magatom_i += 1
                magelem = None
                if self.onlymagpu:
                    magelem = [x for x in self.magdict[mid]["magatom"] if x is not None]
                for dir_ in [FMdir, AFM1dir]:
                    if not os.path.exists(dir_):
                        os.makedirs(dir_)
                    incarfile = os.path.join(dir_, "INCAR")
                    os.system(f"cp {INCARbasic} {incarfile}")
                    if self.pU:
                        INCARPu(element, incar=incarfile, magatomlist=magelem, u_setting=self.puset["u_setting"],
                                openmixparam=self.puset["openmixparam"], mixparam=self.puset["mixparam"])
                    os.system(f"cp {os.path.join(calpath, 'POSCAR' + scalmat)} {os.path.join(dir_, 'POSCAR')}")
                    if os.path.exists(os.path.join(calpath, "POTCAR")):
                        os.system(f"ln -nfs ../POTCAR {os.path.join(dir_, 'POTCAR')}")
                    else:
                        print(f"The POTCAR file of {calpath} does not exist.")

                    with open(incarfile, "a") as f:
                        f.write("MAGMOM = ")
                        if dir_ == FMdir:
                            f.write(" ".join(magmom["FM"]))
                        elif dir_ == AFM1dir:
                            f.write(" ".join(magmom["AFM1"]))
                        f.write("\n")
                if vaspkitKpoints:
                    os.system(
                        f"echo -e \"102\n{kmseshScheme}\n{kmeshrv}\" | (cd {FMdir} && vaspkit >/dev/null 2>&1)")
                    os.system(f"ln -nfs ../FM/KPOINTS {os.path.join(AFM1dir, 'KPOINTS')}")
                self.subgenInputfile(pathlist, ["FM", "AFM1"], scalmat)

            else:
                print(f"Unsupported magnetic state: {magstate}")
                print(f"Please generate the input files manually for {mid}.")

    def catfile(self):
        if len(self.bicalmaterial) != 0:
            with open(os.path.join(self.workdir, "bimaglist.txt"), "w", newline="\n") as f:
                for bimagpath in self.bicalmaterial:
                    f.write(bimagpath + "\n")
            return self.bimaglist
        else:
            print("No magnetic material needs to be calculated again.")
            return ""

    def getEnergy(self, oszipath=None, FMm4=False, copyCONTCAR=True, skipNFM=True):
        """
        Get the energy of the magnetic states.
        After determining the magnetic ground state, the command can be added to add the initial setting of MAGMOM
        in INCAR when performing other task calculations:
        ========================================
        if [ -f MAGMOM ]; then
            cat MAGMOM >> $task/INCAR
        fi
        ========================================

        oszipath: str
            The path of the OSZICAR file.
            If the saveall=True calculation is performed, the OSZICAR file will be saved in a folder named after the
            task, in which case you must specify the task name, e.g. 'opt' or 'scf'.
            If oszipath is "", meaning that the saveall=False calculation is performed, the OSZICAR file will be saved
            in the current folder.
            If oszipath is None, the FM state is skipped.
        FMm4: bool
            If True, the energy of the FM state will be multiplied by 4.
        copyCONTCAR: bool
            Whether to copy the CONTCAR file to the material folder.
        skipFM: bool
            Whether to skip the magstate is "00" (NM/NM) and "03" (NM/FM).
        """
        magdictsave = {}
        bistablemagstate = {}

        if oszipath is None:
            oszipath = self.magpath

        for mid in self.biInputdict.keys():
            magE = {"NFM": [], "AFM1": [], "AFM2": [], "AFM3": []}
            mag = {"NFM": [], "AFM1": [], "AFM2": [], "AFM3": []}
            NFMstatelist = []
            stackmodelist = []
            magstatelist = []
            maglist = []
            magstatelistcorrect = []
            pathlist = self.biInputdict[mid]["submaterial"]
            magstate = self.biInputdict[mid]["magstate"]
            if magstate != "00" and magstate != "03":
                if mid not in magdictsave.keys():
                    magdictsave[mid] = {"stackmode": [], "magstatelist": [], "magstatelistcorrect": [], "maglist": [],
                                        "magE": {"NFM": [], "AFM1": [], "AFM2": [], "AFM3": []},
                                        "mag": {"NFM": [], "AFM1": [], "AFM2": [], "AFM3": []},
                                        "NFMstate": []}
                    # magdictsave[mid]["magstate"] = magstate

                for material in pathlist:
                    magstatestable = ""
                    magEstable = 0
                    magstable = 0
                    outcarpath = None

                    calpath = os.path.join(self.workdir, material)
                    assert os.path.exists(calpath), f"The path {calpath} does not exist."
                    stackmodelist.append("/".join(material.split("/")[1:]))

                    FMdir = os.path.join(calpath, "FM")
                    AFM1dir = os.path.join(calpath, "AFM1")
                    AFM2dir = os.path.join(calpath, "AFM2")
                    AFM3dir = os.path.join(calpath, "AFM3")

                    incarpath = None
                    mayFMm4 = False
                    if os.path.exists(FMdir) and os.path.exists(os.path.join(FMdir, "OSZICAR")):
                        incarpath = os.path.join(FMdir, "INCAR")
                        oszicar = os.path.join(FMdir, "OSZICAR")
                        contcar = os.path.join(FMdir, "CONTCAR")
                    else:
                        mayFMm4 = True
                        if oszipath is not None:
                            oszicar = os.path.join(calpath, oszipath, "OSZICAR")
                            contcar = os.path.join(calpath, oszipath, "CONTCAR")
                        else:
                            oszicar = os.path.join(calpath, "OSZICAR")
                            contcar = os.path.join(calpath, "CONTCAR")
                    if os.path.exists(oszicar):
                        magEfloat = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $5}'",
                                                         shell=True, capture_output=True, text=True).stdout.strip())
                        if FMm4 and mayFMm4:
                            magEfloat *= 4
                        Magfloat = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $10}'",
                                                        shell=True, capture_output=True, text=True).stdout.strip())
                        if abs(Magfloat) <= self.maglowerl:
                            nfmstate_ = "NM"
                        else:
                            nfmstate_ = "FM"
                        NFMstatelist.append(nfmstate_)
                        magE["NFM"].append(magEfloat)
                        mag["NFM"].append(Magfloat)
                        if magEfloat < magEstable:
                            magEstable = magEfloat
                            magstable = Magfloat
                            magstatestable = nfmstate_
                            outcarpath = os.path.join(os.path.dirname(oszicar), "OUTCAR")

                    dirlist = ["AFM1", "AFM2", "AFM3"]
                    for i, dir_ in enumerate([AFM1dir, AFM2dir, AFM3dir]):
                        if os.path.exists(dir_):
                            oszicar = os.path.join(dir_, "OSZICAR")
                            if os.path.exists(oszicar):
                                Magfloat = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $10}'",
                                                                shell=True, capture_output=True,
                                                                text=True).stdout.strip())
                                grepE = subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $5}'",
                                                       shell=True, capture_output=True, text=True)
                                magEfloat = float(grepE.stdout.strip())
                                magE[dirlist[i]].append(magEfloat)
                                mag[dirlist[i]].append(Magfloat)
                                if magEfloat < magEstable:
                                    magEstable = magEfloat
                                    magstable = Magfloat
                                    magstatestable = dirlist[i]
                                    incarpath = os.path.join(dir_, "INCAR")
                                    contcar = os.path.join(dir_, "CONTCAR")
                                    outcarpath = os.path.join(os.path.dirname(oszicar), "OUTCAR")
                    magstatelist.append(magstatestable)
                    maglist.append(magstable)

                    if outcarpath and os.path.exists(outcarpath):
                        if self.Outcar:
                            magmom = self.Outcar(outcarpath).magnetization
                            magmom_list = [magtot['tot'] for magtot in magmom]
                        else:
                            startnu = int(subprocess.run(f"grep -n '# of ion' {outcarpath} | tail -n 1 | cut -d':' -f1",
                                                         shell=True, capture_output=True, text=True).stdout.strip())
                            endnu = int(subprocess.run(f"grep -n 'tot ' {outcarpath} | tail -n 1 | cut -d':' -f1",
                                                       shell=True, capture_output=True, text=True).stdout.strip())
                            nbegin = startnu + 2
                            nend = endnu - 1
                            magmom_list = []
                            for i in range(nbegin, nend):
                                grepmagmom = subprocess.run(f"sed -n '{i}p' {outcarpath}",
                                                            shell=True, capture_output=True, text=True).stdout.strip()
                                magmom_list.append(float(grepmagmom.split()[-1]))
                        magcor = determine_mag_state(magmom_list, self.maglowerl)
                    else:
                        magcor = None
                    magstatelistcorrect.append(magcor)
                    # if magcor is None or magcor == "NM":
                    #     if os.path.exists(os.path.join(calpath, 'MAGMOM')):
                    #         os.system(f"rm {os.path.join(calpath, 'MAGMOM')}")
                    #     incarpath = None

                    if incarpath:
                        grepMAG = subprocess.run(f"grep 'MAGMOM' {incarpath}",
                                                 shell=True, capture_output=True, text=True)
                        bistablemagstate[material] = grepMAG.stdout.strip()
                        if bistablemagstate[material]:
                            with open(os.path.join(calpath, 'MAGMOM'), "w", newline="\n") as f:
                                f.write(f"{bistablemagstate[material]}")

                    if copyCONTCAR and contcar:
                        os.system(f"cp {contcar} {os.path.join(calpath, 'POSCAR')}")

                magdictsave[mid]["stackmode"].extend(stackmodelist)
                magdictsave[mid]["magstatelist"].extend(magstatelist)
                magdictsave[mid]["magstatelistcorrect"].extend(magstatelistcorrect)
                magdictsave[mid]["maglist"].extend(maglist)
                for key in magE:
                    magdictsave[mid]["magE"][key].extend(magE[key])
                for key in mag:
                    magdictsave[mid]["mag"][key].extend(mag[key])
                magdictsave[mid]["NFMstate"].extend(NFMstatelist)

            elif not skipNFM:
                if mid not in magdictsave:
                    magdictsave[mid] = {"stackmode": [], "magstatelist": [], "magstatelistcorrect": [], "maglist": [],
                                        "magE": {"NFM": [], "AFM1": [], "AFM2": [], "AFM3": []},
                                        "mag": {"NFM": [], "AFM1": [], "AFM2": [], "AFM3": []},
                                        "NFMstate": []}
                    # magdictsave[mid]["magstate"] = magstate
                if magstate == "00":
                    stateprefix = "NM"
                else:
                    stateprefix = "FM"
                for material in pathlist:
                    magstatelist.append(stateprefix)
                    calpath = os.path.join(self.workdir, material)
                    assert os.path.exists(calpath), f"The path {calpath} does not exist."
                    stackmodelist.append("/".join(material.split("/")[1:]))

                    mayFMm4 = False
                    FMdir = os.path.join(calpath, "FM")
                    if os.path.exists(FMdir) and os.path.exists(os.path.join(FMdir, "OSZICAR")):
                        oszicar = os.path.join(FMdir, "OSZICAR")
                    else:
                        mayFMm4 = True
                        if oszipath is not None:
                            oszicar = os.path.join(calpath, oszipath, "OSZICAR")
                        else:
                            oszicar = os.path.join(calpath, "OSZICAR")

                    if os.path.exists(oszicar):
                        magEfloat = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $5}'",
                                                         shell=True, capture_output=True, text=True).stdout.strip())
                        if FMm4 and mayFMm4:
                            magEfloat *= 4
                        Magfloat = float(subprocess.run(f"tail -1 {oszicar}" + " | awk '{print $10}'",
                                                        shell=True, capture_output=True, text=True).stdout.strip())
                        if abs(Magfloat) <= self.maglowerl:
                            nfmstate_ = "NM"
                        else:
                            nfmstate_ = "FM"
                        NFMstatelist.append(nfmstate_)
                        magE["NFM"].append(magEfloat)
                        mag["NFM"].append(Magfloat)
                        outcarpath = os.path.join(os.path.dirname(oszicar), "OUTCAR")

                        if os.path.exists(outcarpath):
                            if self.Outcar:
                                magmom = self.Outcar(outcarpath).magnetization
                                magmom_list = [magtot['tot'] for magtot in magmom]
                            else:
                                startnu = int(
                                    subprocess.run(f"grep -n '# of ion' {outcarpath} | tail -n 1 | cut -d':' -f1",
                                                   shell=True, capture_output=True, text=True).stdout.strip())
                                endnu = int(subprocess.run(f"grep -n 'tot ' {outcarpath} | tail -n 1 | cut -d':' -f1",
                                                           shell=True, capture_output=True, text=True).stdout.strip())
                                nbegin = startnu + 2
                                nend = endnu - 1
                                magmom_list = []
                                for i in range(nbegin, nend):
                                    grepmagmom = subprocess.run(f"sed -n '{i}p' {outcarpath}",
                                                                shell=True, capture_output=True,
                                                                text=True).stdout.strip()
                                    magmom_list.append(float(grepmagmom.split()[-1]))
                            magcor = determine_mag_state(magmom_list, self.maglowerl)
                            if magcor == "AFM":
                                magcor = nfmstate_
                            magstatelistcorrect.append(magcor)
                        else:
                            magstatelistcorrect.append(None)

                magdictsave[mid]["stackmode"].extend(stackmodelist)
                magdictsave[mid]["magstatelist"].extend(magstatelist)
                magdictsave[mid]["magstatelistcorrect"].extend(magstatelistcorrect)
                magdictsave[mid]["maglist"].extend(mag["NFM"])
                for mkey, mvalue in magE.items():
                    magdictsave[mid]["magE"][mkey].extend(mvalue)
                for mkey, mvalue in mag.items():
                    magdictsave[mid]["mag"][mkey].extend(mvalue)
                magdictsave[mid]["NFMstate"].extend(NFMstatelist)

        with open(os.path.join(self.workdir, "bimagdict.json"), "w") as f:
            json.dump(magdictsave, f)
        bistablemagstate["comment"] = "if [ ! -f MAGMOM ]; then\n    cat MAGMOM >> $task/INCAR\nfi"""
        with open(os.path.join(self.workdir, "bimagstatecom.json"), "w") as f:
            json.dump(bistablemagstate, f, indent=4)
        summagstate = self.SummagstateDict(magdictsave)
        with open(os.path.join(self.workdir, "bimagdictSum.json"), "w") as f:
            json.dump(summagstate, f, indent=4)

    def SummagstateDict(self, bimagdict):
        """
        Read the dictionary of the magnetic states.
        """
        bimagdictsum = {}
        for mid in bimagdict.keys():
            midprefix = mid.split("_")[0]
            if midprefix not in bimagdictsum.keys():
                bimagdictsum[midprefix] = {"stackmode": [], "magstatelistcorrect": [], "maglist": []}
                # bimagdictsum[midprefix]["magstatelist"] = []
            bimagdictsum[midprefix]["stackmode"].extend(bimagdict[mid]["stackmode"])
            # bimagdictsum[midprefix]["magstatelist"].extend(bimagdict[mid]["magstatelist"])
            bimagdictsum[midprefix]["magstatelistcorrect"].extend(bimagdict[mid]["magstatelistcorrect"])
            bimagdictsum[midprefix]["maglist"].extend(bimagdict[mid]["maglist"])

        for mid in bimagdictsum.keys():
            # stateset = set(bimagdictsum[mid]["magstatelist"])
            # if len(stateset) <= 1:
            #     bimagdictsum[mid]["mixstate"] = False
            # else:
            #     bimagdictsum[mid]["mixstate"] = True
            statesetcorrect = set(bimagdictsum[mid]["magstatelistcorrect"])
            if len(statesetcorrect) <= 1:
                bimagdictsum[mid]["mixstatecorrect"] = False
            else:
                bimagdictsum[mid]["mixstatecorrect"] = True
        return bimagdictsum

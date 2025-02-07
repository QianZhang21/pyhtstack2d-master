import os
from pymatgen.core import Structure


def ciftopos(cifpath=None):
    if cifpath == "" or cifpath is None:
        cifpath = os.getcwd()
        cifpath_basename = ""
    elif not os.path.exists(cifpath):
        print("The path does not exist.")
        raise FileNotFoundError
    else:
        cifpath_basename = os.path.basename(cifpath)
    poscarpath = os.path.join(cifpath_basename, "poscar_dir")
    if not os.path.exists(poscarpath):
        os.makedirs(poscarpath)
    for filename in os.listdir(cifpath):
        if filename.endswith(".cif"):
            cif_path = os.path.join(cifpath, filename)
            structure = Structure.from_file(cif_path)
            poscarfilepath = os.path.join(poscarpath, f"{filename[:-4]}-POSCAR")
            structure.to(fmt="poscar", filename=poscarfilepath)
            print(poscarfilepath)

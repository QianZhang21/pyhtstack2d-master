from pymatgen.core import Structure
import numpy as np
import os
import json

with open("M_Element.json", "r") as f:
    M_dict = json.load(f)

with open("mid.txt", "r") as f:
    lines = f.readlines()

E_BN = -17.754751
save_dict = {}
for li in lines:
    filepath = li.strip().split(" ")[0]
    mag = li.strip().split(" ")[1]
    E0 = li.strip().split(" ")[2]
    stfile = os.path.join(filepath.strip(), "AA", "cord1", "CONTCAR")
    structure = Structure.from_file(stfile)
    a = structure.lattice.a
    b = structure.lattice.b
    gamma = structure.lattice.gamma
    # print(a, b, gamma)
    area = a * b * np.sin(np.radians(gamma))
    labels = structure.labels
    num_BN = len([label for label in labels if label == "B"])
    num_M = len([label for label in labels if label != "B" and label != "N"])
    M_label = [label for label in labels if label != "B" and label != "N"][0]
    E_M = M_dict[M_label]["E0"]
    Ef = float(E0) - E_BN * num_BN - E_M * num_M
    Ef_per_area = Ef / area
    Ef_per_atom = Ef / (num_BN * 2 + num_M)
    Ef_per_aa = Ef_per_area / (num_BN * 2 + num_M)
    # print(f"{filepath.strip()} Mag: {mag}, Ef: {E_f: .4f}, num_BN: {num_BN}, num_M: {num_M}")
    print(f"{filepath.strip()} Mag: {mag}, Ef: {Ef: .4f}, Ef_per_area: {Ef_per_area: .4f}, "
          f"Ef_per_atom: {Ef_per_atom: .4f}, Ef_per_aa: {Ef_per_aa: .6f}")
    save_dict[filepath.strip()] = {"Mag": mag, "Ef": Ef, "Ef_per_area": Ef_per_area,
                                   "Ef_per_atom": Ef_per_atom, "Ef_per_aa": Ef_per_aa,
                                   "num_BN": num_BN, "num_M": num_M, "M_label": M_label,
                                   "a": a, "b": b, "gamma": gamma, "area": area, "E0": E0, "E_M": E_M}

import pandas as pd

df = pd.DataFrame.from_dict(save_dict, orient='index')
df.to_csv("Ef.csv", index=True)

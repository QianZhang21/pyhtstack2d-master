import os

from pymatgen.core.structure import Structure


def build_supercell(stfile, scaling_matrix, save_path=None):
    """
    Builds a supercell from an input structure and a scaling matrix.

    stfile: str
        The path to the input structure file.
    scaling_matrix: list
        The scaling matrix to build the supercell.
    save_path: str
        The path to save the supercell file.
    """
    scalmat = "".join([f"{x}" for x in scaling_matrix])
    if save_path is None:
        save_path = os.path.join(os.path.dirname(stfile), f"{os.path.basename(stfile)}{scalmat}")
    st = Structure.from_file(stfile)
    st.make_supercell(scaling_matrix)
    st.to(filename=save_path, fmt="poscar")
    # try:
    #     st.to(save_path, "poscar")
    # except:
    #     st.to("poscar", save_path)
    # return st


if __name__ == "__main__":
    # build_supercell("../Example2/POSCAR_Gd/1MoS2-POSCAR", [2, 2, 1])
    build_supercell("../tmp/POSCAR_dir0/1NiI2-2.vasp", [2, 2, 1])







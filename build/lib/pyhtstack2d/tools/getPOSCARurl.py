import json
import numpy as np
import os

element_ref_file = os.path.join(os.path.dirname(__file__), 'element_ref.json')
with open(element_ref_file, 'r', encoding='utf-8') as f:
    element_ref_dict = json.load(f)
    element_ref_dict = {int(k): v for k, v in element_ref_dict.items()}

path = os.getcwd()
# if not os.path.isdir(os.path.join(path, "POSCAR_dir")):
#     os.mkdir(os.path.join(path, "POSCAR_dir"))


def WritePOSCAR(data, w_type='C', overwrite=True, lc=None, position_type='C'):
    """
    Get the POSCAR file from the C2DB C1DB, or BIDB.

    data: json
        the json data from the database
    w_type: 'C' or 'F'/'D'
        the type of the POSCAR, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates
    overwrite: bool
        whether to overwrite the existing POSCAR file
    lc: float or int
        the length in the z-direction in a Cartesian coordinate system
        Example: 15.0 la constant is [[...], [...], [0.0, 0.0, 15.0]]
        If "none", the total la constant will be obtained from the database
    position_type: 'C' or 'F'/'D'
        the type of the atomic position in Database, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates
    """
    assert isinstance(data, dict), "The data should be a dictionary."
    assert w_type.upper() in ['C', 'F', 'D'], "The type of the POSCAR is not correct."
    assert position_type.upper() in ['C', 'F', 'D'], "The type of the atomic position is not correct."
    assert lc is None or isinstance(lc, (float, int)), "The length in the z-direction should be a float or an integer."

    w_flag = 'Cartesian' if w_type.upper() == 'C' else 'Direct'
    lattice = data["1"]["cell"]["array"]["__ndarray__"][2]
    atomic_position = data["1"]["positions"]["__ndarray__"][2]

    L_array = np.array(lattice).reshape((3, 3))
    atomic_position = np.array(atomic_position).reshape((-1, 3))
    if position_type.upper() == 'F' or position_type.upper() == 'D':
        P_array = atomic_position.deepcopy()
        atomic_position = np.dot(P_array, L_array)

    if lc is not None:
        assert lc > 0, "The length in the z-direction should be greater than 0."
        L_array[-1, -1] = lc

    if w_type.upper() == 'F' or w_type.upper() == 'D':
        atomic_position = atomic_position @ np.linalg.inv(L_array)

    formula_pos_dict = dict()
    formula_cal_n = data["1"]["numbers"]["__ndarray__"][2]
    n_atomic = data["1"]["numbers"]["__ndarray__"][0][0]
    for i in range(n_atomic):
        atomic_element = element_ref_dict[formula_cal_n[i]]
        if atomic_element not in formula_pos_dict.keys():
            formula_pos_dict[atomic_element] = dict()
            formula_pos_dict[atomic_element]['count'] = 0
            formula_pos_dict[atomic_element]['position'] = []
        formula_pos_dict[atomic_element]['count'] += 1
        formula_pos_dict[atomic_element]['position'].append(atomic_position[i])

    formula_w = ""
    for elem_ in formula_pos_dict.keys():
        formula_w += (elem_ + str(formula_pos_dict[elem_]['count']))

    filename_ = formula_w + "-POSCAR"
    if not os.path.isdir(os.path.join(path, "POSCAR_dir")):
        os.mkdir(os.path.join(path, "POSCAR_dir"))
    save_filename = os.path.join(path, "POSCAR_dir", filename_)
    if not overwrite:
        while os.path.exists(save_filename):
            filename_ = "new_" + filename_
            save_filename = os.path.join(path, "POSCAR_dir", filename_)

    with open(save_filename, "w", encoding='utf-8') as f:
        f.write(formula_w + "\n")
        f.write("1.0\n")
        for i in range(3):
            f.write("    ".join([str(j) for j in L_array[i]]) + "\n")
        for k, v in formula_pos_dict.items():
            f.write(k + "    ")
        f.write("\n")
        for k, v in formula_pos_dict.items():
            f.write(str(v['count']) + "    ")
        f.write("\n")
        f.write(w_flag + "\n")
        for k, v in formula_pos_dict.items():
            for i in v['position']:
                f.write("    ".join([str(j) for j in i]) + "\n")


def DBuidPOSCAR(uid, w_type='C', db="c2db", overwrite=True, lc=None, print_url=False):
    """
    According to the uid of the material in C2DB C1DB or BIDB, get the POSCAR file.

    uid: str or list
        the unique ID of the material in the database
    w_type: 'C' or 'F'/'D'
        the type of the POSCAR, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates
    db: 'C2DB' 'C1DB' or BIDB
        Select the C1DB C2DB or BIDB database
    overwrite: bool
        If true, it overwrites identically named POSCAR files.
    lc: float or int
        the length in the z-direction in a Cartesian coordinate system
        Example: 15.0 la constant is [[...], [...], [0.0, 0.0, 15.0]]
        If "none", the total la constant will be obtained from the database
    print_url: bool
        If true, print the URL of the POSCAR file.
    """
    import json
    import urllib.request

    assert isinstance(uid, (str, list)), "The uid of material should be a string or a list."
    assert w_type.upper() in ['C', 'F', 'D'], "The type of the POSCAR is not correct."
    assert lc is None or isinstance(lc, (float, int)), "The length in the z-direction should be a float or an integer."
    assert db.lower() in ['c2db', 'c1db', 'bidb'] or db is None, "The database should be selected from either C2DB, C1DB or BIDB."

    if isinstance(uid, list):
        for uid_ in uid:
            return DBuidPOSCAR(uid_, w_type, db, overwrite, lc, print_url)

    if "/" in uid:
        uid = uid.split("/")[-1]

    if db.lower() == "c1db":
        prefix_url = 'c1db'
    elif db.lower() == "bidb":
        prefix_url = 'bidb'
    else:
        prefix_url = 'c2db'
    json_position = "https://" + prefix_url + ".fysik.dtu.dk/material/" + uid + '/download/json'

    if print_url:
        print(json_position)

    try:
        with urllib.request.urlopen(json_position) as f:
            data = json.loads(f.read().decode('utf-8'))
            WritePOSCAR(data, w_type, overwrite, lc)
            print(f"The {uid.split('/')[-1]} POSCAR file has been saved in the folder 'POSCAR_dir'.")
    except Exception as e:
        print(f"An error occurred when getting the {uid.split('/')[-1]} POSCAR file.")
        print(f"Please check if the URL '{json_position}' exists. If it does, please try again.\n"
              f"If it does not exist, the uid '{uid.split('/')[-1]}' may not belong to the database.")
        print(e)


def POSCARshiftZ(file_path, overwrite=False, save_type=None, shift=0.5):
    """
    Shift the atomic position of the POSCAR file in z-direction.

    file_path: str
        The path of the POSCAR file
    Overwrite: bool
        If True, the original file will be overwritten.
    save_type: 'C' or None
        The type of the POSCAR, 'C'：Cartesian coordinates; None: Fractional coordinates
    shift: float [0, 1]
        The shift value in the z-direction
    """
    assert os.path.exists(file_path), "The file does not exist."
    assert 0 <= shift <= 1, "The shift value should be in the range [0, 1]."

    with open(file_path, "r", encoding='utf-8') as f:
        data = f.readlines()

    lattice = np.array([i.split() for i in data[2:5]], dtype=float)
    atomic_position = np.array([i.split() for i in data[8:]], dtype=float)
    poscar_type = data[7].lower()

    if poscar_type.startswith('c'):
        atomic_position = atomic_position @ np.linalg.inv(lattice)

    atomic_position[:, -1] += shift
    atomic_position[:, -1] = atomic_position[:, -1] % 1

    if save_type == 'C':
        atomic_position = np.dot(atomic_position, lattice)
        poscar_type = 'Cartesian'
    else:
        poscar_type = 'Direct'

    if overwrite:
        save_filename = file_path
        print(f"The file {file_path} has been overwritten.")
    else:
        save_filename = file_path + "-shift"
        print(f"The file {file_path} has been shifted and saved as {save_filename}.")

    with open(save_filename, "w", encoding='utf-8') as f:
        f.write(data[0])
        f.write(data[1])
        for i in lattice:
            f.write("    ".join([str(j) for j in i]) + "\n")
        f.write(data[5])
        f.write(data[6])
        f.write(poscar_type + "\n")
        for i in atomic_position:
            f.write("    ".join([str(j) for j in i]) + "\n")


import numpy as np
import os


def poscarmkdir():
    """
    Create a directory named 'POSCAR_dir' in the current working directory.
    The directory is used to store the POSCAR files.
    """
    path = os.getcwd()
    if not os.path.isdir(os.path.join(path, "POSCAR_dir")):
        os.mkdir(os.path.join(path, "POSCAR_dir"))
    return os.path.join(path, "POSCAR_dir")


def KeyShow(database, have_gap=True):
    """
    Print the keys of the database.

    have_gap: bool, optional
        Whether the database contains a 'gap' key.
    """
    if have_gap:
        rowone = database.select("gap>0").__next__()
    else:
        rowone = database.get(1)
    key_types = {}
    # print(dir(rowone))
    # print(rowone.get('data'))

    for key in rowone._keys:
        key_type = type(rowone[key]).__name__
        key_types.setdefault(key_type, []).append(key)

    key_types.setdefault("str", []).append("formula")
    key_types.setdefault("int", []).append("natoms")
    for key in ["smax", "mass", "fmax", "volume", "charge"]:
        key_types.setdefault("float", []).append(key)

    for key_type_, keys_ in key_types.items():
        print(f'\"{key_type_}\" keys: \n {keys_}')


def connectDB(database):
    """
    Connect to the C2DB or C1DB database.
    Details of C2DB can be found at: https://cmr.fysik.dtu.dk/c2db/c2db.html
    Details of C1DB can be found at: https://cmr.fysik.dtu.dk/c1db/c1db.html#c1db

    database: str
        Path to the database file.
    """
    from ase.db import connect

    try:
        db_ = connect(database)
        print("Connected to DB.")
        return db_

    except:
        print("Error: Could not connect to C2DB. Please check the database file.", '\n',
              "Details of C2DB can be found at: https://cmr.fysik.dtu.dk/c2db/c2db.html.\n "
              "Details of C1DB can be found at: https://cmr.fysik.dtu.dk/c1db/c1db.html#c1db")
        return None


def WritePOSCARdb(data, overwrite=True, w_type='C', position_type='C', save_path=None):
    """
    Write the POSCAR file from the data.

    data: pos_dict
        The information to write in the POSCAR file.
    overwrite: bool, optional
        Whether to overwrite the existing POSCAR file.
    w_type: 'C' or 'F'/'D'
        the type of the POSCAR, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates.
    position_type: 'C' or 'F'/'D'
        the type of the position in the data, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates.
    save_path: str
        The path to save the POSCAR file.
    """
    assert isinstance(data, dict), "Error: Invalid data type."
    assert w_type.upper() == 'C' or w_type.upper() == 'F' or w_type.upper() == 'D', "Error: Invalid w_type."
    assert position_type.upper() == 'C' or position_type.upper() == 'F' or position_type.upper() == 'D', \
        "Error: Invalid position_type."
    if save_path is not None:
        if not os.path.isdir(save_path):
            os.mkdir(save_path)

    lattice = data['cell']
    positions = data['positions']
    symbols = data['symbols']
    # numbers = data['numbers']
    natoms = data['natoms']
    if natoms is None:
        natoms = len(symbols)
    formula = data['formula']
    if "-" in formula:
        formula = formula.split("-")[0]

    w_flag = 'Cartesian' if w_type.upper() == 'C' else 'Direct'

    if position_type.upper() == 'F' or position_type.upper() == 'D':
        P_array = positions.deepcopy()
        positions = np.dot(P_array, lattice)

    if w_type.upper() == 'F' or w_type.upper() == 'D':
        positions = positions @ np.linalg.inv(lattice)

    formula_pos_dict = dict()
    for i in range(natoms):
        atomic_element = symbols[i]
        if atomic_element not in formula_pos_dict.keys():
            formula_pos_dict[atomic_element] = dict()
            formula_pos_dict[atomic_element]['count'] = 0
            formula_pos_dict[atomic_element]['position'] = []
        formula_pos_dict[atomic_element]['count'] += 1
        formula_pos_dict[atomic_element]['position'].append(positions[i])

    filename_ = formula + "-POSCAR"
    save_filename = os.path.join(save_path, filename_)
    if not overwrite:
        while os.path.exists(save_filename):
            filename_ = "new_" + filename_
            save_filename = os.path.join(save_path, filename_)
    with open(save_filename, "w", encoding='utf-8', newline='\n') as f:
        f.write(formula + '\n')
        f.write('1.0\n')
        for i in range(3):
            f.write('    '.join([str(j) for j in lattice[i]]) + '\n')

        for key in formula_pos_dict.keys():
            f.write(key + '    ')
        f.write('\n')
        for key in formula_pos_dict.keys():
            f.write(str(formula_pos_dict[key]['count']) + '    ')
        f.write('\n')
        f.write(w_flag + '\n')
        for key in formula_pos_dict.keys():
            for pos in formula_pos_dict[key]['position']:
                f.write('    '.join([str(i) for i in pos]) + '\n')


def getuid(database, criteria, have_gap=True, save_csv=False, save_csv_path=None):
    """
    Retrieve specific rows from C2DB or C1DB based on selection criteria.

    database: str
        Path to the database file.
    criteria: str
        Selection criteria for the database. For example, 'gap>0, stoichiometry=ABC'.
    key_show: bool, optional
        Whether to print the keys of the database.
    have_gap: bool, optional
        Whether the database contains a 'gap' key, only the key_show is True, this parameter is useful.
    save_csv: bool, optional
        Whether to save the retrieved data in a CSV file.
    save_csv_path: str, optional
        The path to save the CSV file.

    """
    db_ = connectDB(database)
    try:
        rows = db_.select(criteria)
        uid_list = []
        if save_csv:
            import pandas as pd
            key_value_pairs_list = []
            for row in rows:
                key_value_pairs_list.append(row.key_value_pairs)
                uid_list.append(row.get("uid"))
            key_value_pairs_df = pd.DataFrame(key_value_pairs_list)

            key_value_pairs_df = key_value_pairs_df.drop(columns=['folder'])
            if save_csv_path:
                try:
                    key_value_pairs_df.to_csv(save_csv_path)
                    print("CSV file saved at:", save_csv_path)
                except:
                    print("Error: Could not save the CSV file at the specified path. Please check the path.")
            else:
                key_value_pairs_df.to_csv("rows.csv")
                print("CSV file saved at: rows.csv")
        else:
            for row in rows:
                uid_list.append(row.uid)
        return uid_list
    except:
        print("Error: Could not retrieve the specific rows from C2DB or C1DB. Please check the selection criteria.")
        print("\nThe selectable keys are:")
        KeyShow(db_, have_gap)
        print("\nThe selection criteria should be in the format: 'gap>0, stoichiometry=ABC'")
        return None


def getPOSACRfromuid(database, uid, overwrite=True, w_type='C', position_type='C', save_path=None):
    """
    Get the POSCAR file from the database based on the uid.

    database: str
        Path to the database file.
    uid: str
        The uid of the data.
    overwrite: bool, optional
        Whether to overwrite the existing POSCAR file.
    w_type: 'C' or 'F'/'D'
        the type of the POSCAR, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates.
    position_type: 'C' or 'F'/'D'
        the type of the position in the data, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates.
    save_path: str, optional
        The path to save the POSCAR file.
    """
    db_ = connectDB(database)
    if not save_path:
        save_path = poscarmkdir()
    else:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
    try:
        if isinstance(uid, str):
            uid = [uid]
        for uid_ in uid:
            row = db_.get(uid=uid_)
            data_dict = {'cell': row.get('cell'), 'positions': row.get('positions'), 'symbols': row.get('symbols'),
                         'natoms': row.get('natoms'), 'formula': row.get('uid')}
            WritePOSCARdb(data_dict, overwrite=overwrite, w_type=w_type, position_type=position_type, save_path=save_path)
        print(f"POSCAR file saved in {save_path}")
    except:
        print("Error: Could not retrieve the specific rows from C2DB or C1DB. Please check the selection criteria.")
        print("\nThe selectable keys are:")
        KeyShow(db_, have_gap=True)
        print("\nThe selection criteria should be in the format: 'gap>0, stoichiometry=ABC'")


def getPOSACR(database, criteria, overwrite=True, w_type='C', position_type='C', save_path=None):
    """
    Retrieve specific rows from C2DB or C1DB based on selection criteria, and save the POSCAR files.

    database: str
        Path to the database file.
    criteria: str
        Selection criteria for the database. For example, 'gap>0, stoichiometry=ABC'.
    overwrite: bool, optional
        Whether to overwrite the existing POSCAR file.
    w_type: 'C' or 'F'/'D'
        the type of the POSCAR, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates.
    position_type: 'C' or 'F'/'D'
        the type of the position in the data, 'C'：Cartesian coordinates; 'F'/'D': Fractional coordinates.
    save_path: str, optional
        The path to save the POSCAR file.
    """
    db_ = connectDB(database)
    if not save_path:
        save_path = poscarmkdir()
    else:
        if not os.path.exists(save_path):
            os.makedirs(save_path)

    try:
        rows = db_.select(criteria)
        for row in rows:
            data_dict = {'cell': row.get('cell'), 'positions': row.get('positions'), 'symbols': row.get('symbols'),
                         'natoms': row.get('natoms'), 'formula': row.get('uid')}
            WritePOSCARdb(data_dict, overwrite=overwrite, w_type=w_type, position_type=position_type, save_path=save_path)
        print(f"POSCAR file saved in {save_path}")
    except:
        print("Error: Could not retrieve the specific rows from C2DB or C1DB. Please check the selection criteria.")
        print("\nThe selectable keys are:")
        KeyShow(db_, have_gap=True)
        print("\nThe selection criteria should be in the format: 'gap>0, stoichiometry=ABC'")


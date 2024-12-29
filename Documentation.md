# 1. What is PyHTStack2D?

**PyHTStack2D** is a Python package developed to support high-throughput computational research on two-dimensional (2D) van der Waals (vdWs) stacking structures. This package addresses the complexities associated with high-throughput stacking layered vdWs structures by providing several key functionalities:

### 1.1. Batch Stacking of Structures

PyHTStack2D facilitates the systematic stacking of both homo- and heterostructures. It incorporates specific stacking sequences and patterns, such as those seen in phase transitions in transition metal dichalcogenides (TMDs), allowing researchers to efficiently generate diverse stacking configurations.

### 1.2. Automated Directory and Script Generation

The package automates the creation of computational directories and batch submission shell scripts, streamlining high-throughput calculations and reducing manual setup time.

Additionally, PyHTStack2D allows users to extract POSCAR files for 2D material structures from the well-known Computational 2D Materials Database (C2DB) based on specific filtering criteria. It also supports batch extraction and preliminary analysis of results after high-throughput computations, providing an efficient way to explore and analyze 2D stacking materials.

# 2. Configuration Requirements

To successfully execute the intended simulations, the system should adhere to the following configuration prerequisites:

### Software Dependencies:
- **Python**: A minimum of version 3.0 is required to support the necessary computational libraries and scripts.
- **Shell Environment**: The Bash shell is essential for scripting and managing the sequence of commands that initiate and control simulations.

### Computational Software:
The system should be equipped with one or more of the following computational chemistry and physics software packages, essential for conducting electronic structure calculations:
- **VASP (Vienna Ab initio Simulation Package)**: This package requires a license for both academic and commercial use and offers tools for molecular dynamics simulations and electronic structure calculations based on density-functional theory (DFT). [VASP website](https://www.vasp.at)
- **Quantum ESPRESSO (QE)**: An open-source suite of tools designed for electronic-structure calculations using plane wave basis sets. [Quantum ESPRESSO website](https://www.quantum-espresso.org)
- **PWDFT (Plane Wave Density Functional Theory)**: A first-principles calculation software that solves the Kohn-Sham equations using plane wave basis sets, ideal for periodic system analyses. DOI: [10.1021/acs.jctc.6b01184](https://doi.org/10.1021/acs.jctc.6b01184)


# 3. Installation

To install **PyHTStack2D**, you have two options:

### Option 1: Install from Source

1. Download `pyhtstack2d-master.tar.gz` and extract the contents.

2. In the extracted directory, install the package by running:
   
   ```bash
   python setup.py install
   ```

### Option 2: Install from Wheel

1. Download `dist/pyhtstack2d-0.1-py3-none-any.whl` and `dist/install_with_config.py`.

2. Run the following command to install with configuration:
   
   ```bash
   python install_with_config.py
   ```
   
   During installation, you will be prompted to specify the directory path for `POTCAR` files and configure the machine settings for job submission.

> **Note:** If you install the package directly using the wheel file (`pyhtstack2d-0.1-py3-none-any.whl`) without `install_with_config.py`, you may need to manually create a configuration file at `~/.config/.PyHTStack2D.json`.

After installation, you can manually edit parameters in `~/.config/.PyHTStack2D.json` at any time to adjust settings.

### Requirements

The **PyHTStack2D** package depends on several libraries. Below is a list of required dependencies and their recommended versions:

- `ase==3.22.1`
- `matplotlib==3.5.1`
- `numpy==2.1.2`
- `pandas==1.4.4`
- `pymatgen==2023.11.12`
- `vaspkit>=1.3.5` (Note: **vaspkit** is not a Python package and needs to be downloaded and installed separately.)

During installation, you will be prompted with:

```plaintext
Do you want to install recommended dependencies? (yes/no): 
```

- If you choose `yes`, the installer will automatically install the recommended versions of these dependencies.
- If you choose `no`, you may need to install each dependency manually or use different compatible versions of these libraries as per your system setup.

For `vaspkit`, please follow the official installation instructions on the [Vaspkit website](https://vaspkit.com) to install it separately.


# 4. Modules

## 4.1 POSCAR retrieval from C2DB and processing

### 4.1.1. `pyhtstack2d.tools.queryDB`

This module provides functions to query the Computational 2D Materials Database (C2DB) and extract 2D material structures based on specific filtering criteria.

**Use cases**:

```python
# Define specific criteria for filtering 2D materials from C2DB.
criteria = 'natoms=3, dyn_stab=Yes, hform<0, gap>0, layergroup=p3m1'
# Database file obtained from C2DB developers.
db = "c2db.db"

# Retrieve POSCAR files for 2D materials from C2DB based on specified criteria.
from pyhtstack2d.tools.queryDB import getPOSACR
getPOSACR(database=db, criteria=criteria)

# Fetch the unique ID of 2D materials from C2DB based on specific criteria, and save the results to a CSV file.
from pyhtstack2d.tools.queryDB import getuid
getuid(database=db, criteria=criteria, save_csv_path='tmd_p3m1.csv', save_csv=True)

# Filter for direct bandgap 2D materials from the CSV file.
# Keywords 'uid', 'gap', and 'gap_dir' are as defined in C2DB.
import pandas
tmd_p3m1 = pandas.read_csv('tmd_p3m1.csv')
tmd_p3m1_uid_dir = []
for i in range(len(tmd_p3m1)):
    if abs(tmd_p3m1.iloc[i]['gap'] - tmd_p3m1.iloc[i]['gap_dir']) < 1e-3:
        tmd_p3m1_uid_dir.append(tmd_p3m1.iloc[i]['uid'])

# Retrieve POSCAR files for selected 2D materials from C2DB using their unique IDs.
from pyhtstack2d.tools.queryDB import getPOSACRfromuid
for uid_i in tmd_p3m1_uid_dir:
    getPOSACRfromuid(database=db, uid=uid_i)
```

### 4.1.2 `pyhtstack2d.tools.queryDBurl`

This module retrieves the POSCAR file for a specified material directly from a webpage URL using its unique identifier (uid).

**Use cases**:

```python
from pyhtstack2d.tools.getPOSCARurl import DBuidPOSCAR
# Specify the unique ID of the material
uid = '1MoS2-1'
# Retrieve the POSCAR file from C2DB using the uid, with the option to overwrite any existing file
DBuidPOSCAR(uid=uid, db="c2db", overwrite=True)
```

### 4.1.3 `pyhtstack2d.tools.cif2pos`

This module facilitates the conversion of CIF files to POSCAR format.

**Use cases**:

```python
from pyhtstack2d.tools.cif2pos import ciftopos
# Specify the directory path containing the CIF files; set to None if the files are in the current directory.
# This command will convert all *.cif files in the specified directory to POSCAR format.
ciftopos(cifpath='./cif_dir')
```

### 4.1.3 `pyhtstack2d.tools.POSCARelemc`

This module is used to modify the elements in a POSCAR file.

**Use cases**:

```python
from pyhtstack2d.tools.POSCARelemc import modify_elements_in_poscar
# Specify the POSCAR file path and the new elements to replace the existing ones.
file_path = "Ag2ReBr6-POSCAR"
new_elements = ['Ag', 'V', 'P', 'S']
save_path = "POSCAR_dir"
# Modify the elements in the POSCAR file and save the modified file to the specified directory.
modify_elements_in_poscar(file_path, new_elements, save_path)
```

### 4.1.4 `pyhtstack2d.tools.renameposcar`

For stacking purposes, POSCAR filenames need to follow the "uid-POSCAR" format. This module renames POSCAR files that are in the "POSCAR-uid" format.

**Use cases**:

```python
from pyhtstack2d.tools.renameposcar import swap_poscar_filename
# Specify the folder containing POSCAR files. This function renames files as needed to match the "*-POSCAR" format.
swap_poscar_filename(folder_path="POSCAR_dir")
```

### 4.1.5 `pyhtstack2d.buildbilayer.builsupercell`

This module leverages [Python Materials Genomics (Pymatgen)](https://pymatgen.org/) to create a supercell from a specified POSCAR file.

**Use Case**:

```python
from pyhtstack2d.buildbilayer.builsupercell import build_supercell
# Specify the POSCAR file path and the scaling matrix for the supercell.
stfile = "MoS2-POSCAR"
scaling_matrix = [2, 2, 1]
# Generate the 2x2x1  supercell.
build_supercell(stfile=stfile, scaling_matrix=scaling_matrix)
```

## 4.2. Stacking Bilayer Structures

### 4.2.1 `pyhtstack2d.buildbilayer.stackBilayer`

This module is designed to facilitate the stacking of two 2D materials to form a bilayer structure. Currently, seven classes are available for constructing bilayer configurations:

- **`Bilayer`**: A fundamental bilayer stacking class suitable for general 2D material combinations.
- **`TMDHBilayer`**: Constructs bilayer structures for transition metal dichalcogenides, halides, or oxides (TMDs/MHs/TMOs).
- **`MNXYBilayer`**: Enables the construction of bilayer structures for materials with Al₂S₂-like configurations, featuring layer groups P-3m1 (1T) or P-6m2 (2H).
- **`N2Bilayer`**: Creates bilayer structures where each layer is composed of two atoms (N₂ configuration).
- **`N2TMDHbilayer`**: Constructs a bilayer with one layer as an N₂ monolayer and the other as a TMD/MH/TMO monolayer (e.g., Graphene-MoS₂). Currently, this class supports only 2H-type TMDs, MHs, and TMOs.
- **`N2MNXYBilayer`**: Enables construction of bilayers where one layer is an N₂ monolayer, and the other is an MNXY monolayer (e.g., Graphene-Al₂S₂). This class currently supports only 2H-type MNXY materials.
- **`TMDHBilayerSquare`**: Constructs square-lattice TMD/MH/TMO bilayer structures with layer group p-4m2 symmetry.

These specialised classes support a variety of highly symmetric stacking configurations, making them adaptable to a wide range of 2D bilayer structure requirements.

**Parameters:**

- **`st1` or `st2`**: List of elements or POSCAR path of monolayer. When this is a POSCAR path, other structural settings (e.g., `la1`, `position1`) can be ignored.
- **`la1` or `la2`**: Lattice matrix for each monolayer when elements are provided as lists and class is `Bilayer`. Otherwise, this is the float value of the lattice constant.
- **`position1` or `position2`**: Atomic positions of each monolayer when using element lists.
- **`d_intra1` or `d_intra2`**: For non-`Bilayer` classes, this specifies the z-direction intralayer distance.
- **`d_inter`**: Initial interlayer distance between the monolayers.
- **`lv`**: Vacuum layer thickness.
- **`mismatch_threshold`**: Maximum mismatch threshold for the lattice constants of the two monolayers.
- **`savepath`**: Path to save the POSCAR file. If `None`, the POSCAR is saved in the default `BiPOSCAR` directory.
- **`stackmode`**: Stacking mode for the two monolayers, used when class is `Bilayer`. Defaults to "AA".
- **`savenamemode`**: Naming mode for saved POSCAR files.
  - **1**: Files saved in `"{savepath}/{formula}_{mismatch}/{genmode}/{cord*}/POSCAR"`.
  - **2**: Files saved as `"{savepath}/{formula}_{genmode}_{cord*}-{mismatch}-POSCAR"`.
  
**Use Cases:**

For instance, you can stack MoS₂ and WSe₂ bilayer structures by specifying the elements, lattice parameters, and atomic positions, or by using POSCAR files.

```python
from pyhtstack2d.buildbilayer.stackBilayer import Bilayer

# Define elements, lattice, and positions for MoS₂
elem1 = ['Mo', 'S', 'S']
la1 = [[3.184, 0.0, 0.0], [-1.592, 2.757, 0.0], [0.0, 0.0, 18.127]]
pos1 = [[0.0, 0.0, 0.50], [2/3, 1/3, 0.59], [2/3, 1/3, 0.41]]
# Define elements, lattice, and positions for WSe₂
elem2 = ['W', 'Se', 'Se']
la2 = [[3.319, 0.0, 0.0], [-1.659, 2.874, 0.0], [0.0, 0.0, 18.358]]
pos2 = [[0.0, 0.0, 0.50], [2/3, 1/3, 0.59], [2/3, 1/3, 0.41]]
# Construct AA stacking of MoS₂/WSe₂
Bilayer(elem1, elem2, la1=la1, la2=la2, position1=pos1, position2=pos2).WritePOSCAR()

# Construct AB stacking of MoS₂/WSe₂
pos2 = [[2/3, 1/3, 0.50], [0.0, 0.0, 0.59], [0.0, 0.0, 0.41]]
Bilayer(st1=elem1, st2=elem2, la1=la1, la2=la2, position1=pos1, position2=pos2, stackmode="AB").WritePOSCAR()

# Alternatively, stack using POSCAR files
poscar1 = "MoS2-POSCAR"
poscar2 = "WSe2-POSCAR"
Bilayer(st1=poscar1, st2=poscar2).WritePOSCAR()

# Create a three-layer stack of MoS₂/WSe₂/MoS₂
biposcar = "BiPOSCAR dir/S2MoSe2W_4.15%_3/AA/cord1/POCAR"
Bilayer(st1=biposcar, st2=poscar1, d_inter=3.5, lv=35.0).WritePOSCAR()
```

For non-`Bilayer` classes, structural information should be entered sequentially according to atomic positions in the z-direction, for example:

```python
from pyhtstack2d.buildbilayer.stackBilayer import TMDHBilayer
elem1 = ['S', 'Mo', 'S']
a1 = 3.184
d_intra1 = [1.63, 1.63]
TMDHBilayer(st1=elem1, a1=a1, d_intra1=d_intra1).WritePOSCAR()
```

### 4.2.2. `pyhtstack2d.buildbilayer.batchStackBilayer`

This module provides a class for batch stacking of bilayer structures.

**Parameters:**

- **`pos_dir` or `pos_dir2`**: Directory containing POSCAR files for monolayers. When only `pos_dir` is specified, each POSCAR file within the directory is paired with every other POSCAR in the directory for stacking. When both `pos_dir` and `pos_dir2` are specified, each POSCAR in `pos_dir` will be stacked with each POSCAR in `pos_dir2`.
- **`la_mismatch` or `lb_mismatch`**: Maximum allowable mismatch threshold for the lattice constants `a` and `b`.
- **`homo`**: Determines whether to stack the same monolayer with itself. Defaults to `False`.
- **`genmode`**: Specifies the class to use for stacking. Options include:
  - `"bilayer"`: Uses the `Bilayer` class
  - `"tmdh"`: Uses the `TMDHBilayer` class
  - `"mnxy"`: Uses the `MNXYBilayer` class
  - `"n2"`: Uses the `N2Bilayer` class
  - `"n2tmdh"`: Uses the `N2TMDHBilayer` class
  - `"n2mnxy"`: Uses the `N2MNXYBilayer` class
  - `"tmdhsquare"`: Uses the `TMDHBilayerSquare` class
- **`kwargs`**: Additional parameters for the specified stacking class.

**Use Case:**

```python
from pyhtstack2d.buildbilayer.batchStackBilayer import GenBiLayer
# The following example demonstrates batch stacking with specified mismatch tolerance, mode, and stacking of identical monolayers (homo-structures).
GenBiLayer(pos_dir="POSCAR_dir",la_mismatch=3.0,homo=True,genmode="tmdh").batch_stack()
```


### 4.2.3. `pyhtstack2d.buildbilayer.RotMovePOSCAR`

This module provides two primary functions for rotating and moving atomic positions.

Combining this module with the bilayer stack generation mode enables the creation of stacking results similar to the `tmdh` generation mode.

**Use Case:**

```python
from pyhtstack2d.buildbilayer.RotMovePOSCAR import rotate_poscar, move_poscar

# Path to the POSCAR file of a monolayer MoS2
poscar = "POSCAR/1MoS2-POSCAR"

# Move the atomic positions of the POSCAR file. The moved POSCAR files will be saved in the directory "POSCAR_moved".
move_poscar(poscar, [0, 0, 0])            # No translation; original position
move_poscar(poscar, [1/3, -1/3, 0])       # Translate by (1/3, -1/3, 0) in fractional coordinates
move_poscar(poscar, [-1/3, 1/3, 0])       # Translate by (-1/3, 1/3, 0) in fractional coordinates

# Path to the POSCAR file of a monolayer WS2
poscar2 = "POSCAR/1WS2-POSCAR" # or "POSCAR/2HfIN-2.vasp"

# Rotate the atomic positions of the POSCAR file. The rotated POSCAR files will be saved in the directory "POSCAR_rotated". 
rotate_poscar(poscar2, 0)                 # No rotation; original orientation
rotate_poscar(poscar2, 60)                # Rotate by 60 degrees

from pyhtstack2d.buildbilayer.batchStackBilayer import GenBiLayer

# Generate bilayer stacking configurations using the modified POSCAR files
# overwrite: If False, do not overwrite existing stacked POSCAR files
# skip_xy_rev: If True, skip automatic consideration of monolayer rotation for stacking
GenBiLayer(pos_dir="POSCAR_moved", pos_dir2="POSCAR_rotated", genmode="bilayer",
           overwrite=False, skip_xy_rev=True).batch_stack()
```

> **Note:** The rotation is a pseudo-rotation that only rotates the atomic coordinates and only facilitates stacking. For example, hexagonal cells can only be rotated 60, 120 and 180 degrees, other angles will change the original structure.

### 4.2.4. `pyhtstack2d.buildbilayer.TransLattice`

This module provides a class for transforming the lattice, such as converting a hexagonal structure to an orthorhombic structure, or constructing a non-integer supercell.

#### Parameters:

- **`filename`**: Path to the POSCAR file of the hexagonal structure.
- **`transmat`**: Transformation matrix for converting the cell. Defaults to `[[-1, 1, 0], [-1, -1, 0], [0, 0, 1]]`.
- **`overwrite`**: Determines whether to overwrite the original POSCAR file. Defaults to `False`.

#### Use Case:

```python
from pyhtstack2d.buildbilayer.TransLattice import transformlattice

# For the 2H-MoS2 system is characterised by a hexagonal structure with a 120 degree angle between the lattice vectors a and b.
# Convert the hexagonal structure to an orthorhombic structure with the smallest lattice constants.
transformlattice(filename="MoS2-POSCAR", transmat=None, overwrite=False)

import numpy as np
# Custom transformation matrix, which can be used to convert the hexagonal structure to an orthorhombic structure.
transmat = np.array([[-2, 0, 0], [-1, -2, 0], [0, 0, 1]])
transformlattice(filename="MoS2-POSCAR", transmat=transmat)

# sqrt(7) x sqrt(7) x 1 supercell of 2H-MoS2
transmat = np.array([[3, 2, 0], [-2, 1, 0], [0, 0, 1]])
transformlattice(filename="MoS2-POSCAR", transmat=transmat)
```

### 4.2.5. `pyhtstack2d.buildbilayer.RemvDuplicates`

This module provides a class for removing duplicate structures from a directory of POSCAR files.

#### Parameters:

- **`posdir`**: Directory containing the POSCAR files.
- **`multist`**: Flag indicating whether there are multiple different systems in the directory. Defaults to `True`.
- **`threshold`**: Tolerance for comparing the structures. Defaults to `1e-2`.
- **`deleteempty`**: Flag indicating whether to delete empty directories. Defaults to `True`.

#### Example Usage:

```python
from pyhtstack2d.buildbilayer.RemvDuplicates import remove_duplicates

# Remove duplicate structures from the directory "BiPOSCAR_dir"
remove_duplicates("BiPOSCAR_dir")
```


## 4.3. Generating Input Files

### Fundamental Parameters

- **`stfile`**: Specifies a single POSCAR file name or a list name of POSCAR files.
- **`posdir`**: Directory where the POSCAR files are stored. If `stfile` is not specified, the program will iterate through the folders in `posdir`. The default is `"POSCAR_dir"`.
- **`workdir`**: Directory where the generated input files are stored. The default directory is `"work"`.
- **`user_settings` or `user_incar_settings`**: User-defined INCAR settings for calculations.
- **`overwrite`**: Determines whether to overwrite existing files.
- **`taskname`**: Name used as a suffix or prefix for input filenames in multi-task calculations.

These fundamental parameters are referenced across modules 4.3.1 to 4.3.5.

### 4.3.1. `pyhtstack2d.calcSets.qeSets`

This module generates input files for Quantum Espresso (QE) calculations based on POSCAR files.

**Special Parameters**

- **`spinOpen`**: Indicates if the input file should be configured for spin-polarized calculations.
- **`inputOpt`**: Specifies the calculation template to generate, such as `'scf'`, `'relax'`, or `'md'`.
- **`block_setting`**: Controls the settings for modifying block tags and ending markers in the input file. The default configuration for `block_setting` is as follows: `{"control": {"blockname": "&CONTROL", "end": "/"}, ...}`

**Use Case**

```python
from pyhtstack2d.calcSets.qeSets import InputWriter4qe
# Generating a Quantum Espresso input file with custom settings
InputWriter4qe(workdir="QE_dir", inputOpt="relax", user_settings={"ecutwfc": 50, "ecutrho": 400, "kpoints": [6, 6, 1]}).write_input()
```

### 4.3.2. `pyhtstack2d.calcSets.pwdftSet`

This module generates input files for PWDFT (Plane-Wave Density Functional Theory) calculations based on POSCAR files.

**Special Parameters**

- **`potpath`**: Specifies the directory path to the pseudopotential files required for the calculation.
- **`pseudo_suffix`**: Defines the suffix for the pseudopotential filenames.

**Use Case**

```python
from pyhtstack2d.calcSets.pwdftSet import InputWriter4pwdft
# Generating a PWDFT input file for MoS2 with a specified pseudopotential suffix
InputWriter4pwdft(stfile='MoS2-POSCAR',workdir="PWDFT_dir",pseudo_suffix="_ONCV_PBE-1.0.upf").write_input()
```

### 4.3.3. `pyhtstack2d.calcSets.vaspSetsWriter`

This module generates input files for VASP (Vienna Ab initio Simulation Package) calculations based on POSCAR files.

**Special Parameters**

- **`inputOpt`**: Specifies the calculation mode, with options including `'relax'`, `'scf'`, `'band'`, `'optic'`, and `'eps'` (electrostatic potential). Additional suffixes denote specific functional or correction types:
  - `-hse06`: HSE06 hybrid functional
  - `-pu`: DFT+U
  - `-dip`: Dipole correction
  - `<-d2, -d3, -d3bj, -optpbe, -optb88, -optb86b>`: van der Waals (vdW) corrections
  - For instance, `'scf-d3'` indicates a static calculation with DFT-D3 no-damping correction. If `inputOpt='VePu'`, an INCAR file is generated for DFT+U and valence electron calculations.
- **`ini_magmom`**: Sets the initial magnetic moment, if required.
- **`magmom_setting`**: User-defined settings for initial magnetic moments.
- **`nonmagmom`**: Sets the initial magnetic moment for non-magnetic (NM) atoms, with a default of 0.
- **`kselfSet`**: Custom KPOINTS settings, e.g., `kselfSet=[11, 11, 1]`.
- **`kmeshrv`**: Defines the k-mesh accuracy level.
- **`gamma`**: Determines if the gamma-centered k-point should be used.
- **`is2D`**: Specifies if the material is 2D.
- **`postype`**: Used for generating the KPOINTS file for band structure calculations. If `postype` is specified and `inputOpt` includes `'band'`, a line-mode k-path is generated.
- **`angle`**: Sets the angle of the hexagonal lattice for line-mode k-path generation; values can be 60 or 120.
- **`potpath`**: Path to the pseudopotential file.
- **`is_print`**: If enabled, prints the path where the input files are saved.

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspSetsWriter import InputWriter
import os
from glob import glob
# Example: Generate VASP input files with DFT-D3 no-damping correction for a static calculation
filenames = glob(os.path.join("POSCAR_dir", "*"))
InputWriter(filenames, inputOpt='scf-d3-dip', kmeshrv=0.04).write_input()

# Example: For nested directories and multi-task calculations
# When POSCAR files are stored in nested directories, use `write_input_multi_pos()` instead of `write_input()`
InputWriter(posdir='BiPOSCAR_dir',inputOpt='scf-d3-dip',kmeshrv=0.04,taskname='scf').write_input_multi_pos()
```

### 4.3.4. `pyhtstack2d.calcSets.vaspPmgWriter`

This module generates input files for VASP calculations using the pymatgen interface, based on POSCAR files.

**Special Parameters**

- **`inputOpt`**: Defines the calculation mode, with options corresponding to pymatgen classes. For example, options include `"MPRelax"` and other modes as per the pymatgen documentation.
- **`is_print`**: If enabled, prints the path where the input files are saved.

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspPmgWriter import PmgInputWriter
# Display available setting modes
PmgInputWriter().mode_show()
# Example: Generate input files with a specific pymatgen mode and task name suffix
# Iterates through files in the default "POSCAR_dir" folder
PmgInputWriter(inputOpt="MPStatic", taskname="scf").write_input()

# Example: For nested directories and multi-task calculations
# When POSCAR files are stored in nested directories, use `write_input_multi_pos()` instead of `write_input()`
PmgInputWriter(posdir="BiPOSCAR_dir").write_input_multi_pos()
```

### 4.3.5. `pyhtstack2d.calcSets.vaspVaspkitWriter`

This module generates input files for VASP calculations using the Vaspkit interface, based on POSCAR files. Typically, this module creates a `mkdir.sh` script that needs to be executed to complete the full setup of input files.

**Special Parameters**

- **`inputOpt`**: Specifies the calculation mode, aligning with Vaspkit's INCAR options menu. For example, `inputOpt = 'STH6D3'` represents a static calculation with HSE06-D3.
- **`incexis`, `potexis`, `kpexis`**: Flags to skip the generation of INCAR, POTCAR, or KPOINTS files if they already exist, useful for scenarios where input files are shared across structures.
- **`kmseshScheme`**: Defines the k-mesh scheme:
  - 1: Monkhorst-Pack Scheme
  - 2: Gamma-Centered Scheme
  - 3: Irreducible K-Points with Gamma Scheme
- **`kmeshrv`**: Sets the k-mesh resolved value:
  - Accuracy levels: Gamma-Only (0), Low (0.06–0.04), Medium (0.04–0.03), Fine (0.02–0.01)
- **`bandkpath`**: K-path option for band structure calculations:
  - 301: 1D structure
  - 302: 2D structure
  - 303: 3D structure
- **`Hybridbandkpath`**: Configures the k-path for hybrid band structure calculations. For example, `[251, 2, 0.03, 0.04]` sets up a Gamma-centered scheme with SCF calculation k-mesh resolved at 0.03 and band calculation at 0.04.
- **`genpotcar`**: Specifies the pseudopotential generation scheme:
  - 103: Default POTCAR generation
  - 104: User-specified POTCAR generation
- **`shoverwrite`**: Determines whether to overwrite an existing `mkdir.sh` file.
- **`is_print`**: If enabled, prints the path where input files are saved.

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspVaspkitWriter import VaspkitInputWriter
from glob import glob
import os
# Example: Generate input files for HSE06-D3 static calculation with a hybrid band structure k-path setting
filenames = glob(os.path.join("POSCAR_dir", "*"))
VaspkitInputWriter(filenames,inputOpt="STH6D3",Hybridbandkpath=[251, 1, 0.04, 0.04],taskname="hybband").write_input()

# Example: For nested directories
# Use `write_input_multi_pos()` for POSCAR files in nested directories
VaspkitInputWriter(posdir="BiPOSCAR_dir",taskname="scf").write_input_multi_pos()
```

### 4.3.6. `pyhtstack2d.tools.genInput.GenRunDir`

This module is designed for batch generation of input files for VASP calculations, with flexible settings for handling various submission and directory structures.

**Parameters**

- **`posdir`**: Directory where the POSCAR files are stored. The default is `"POSCAR_dir"`.
- **`posname`**: Specifies character patterns in the POSCAR filenames. For instance, `"*"` selects all files, while `"*POSCAR*"` selects files containing "POSCAR" in their names.
- **`workdir`**: Directory where the generated input files are stored. 
- **`multilevel`**: Indicates the number of sub-directory levels used for writing the `run.sh` file for multiple POSCAR files organized by chemical formula. If `multilevel` is specified, `posname` and `workdir` are ignored. For example, if the POSCAR files are stored in `"BiPOSCAR_dir/S24P8Ag4V4_0.00%/AA/cord1"`, `multilevel` would be `4`.
- **`subset`**: Determines the job submission system:
  - `"pbs"`: PBS submission mode
  - `"slurm"`: SLURM submission mode
- **`mpirun`**: Specifies the MPI command.
- **`vasp`**: Name of the VASP module to use, e.g., `"vasp_std"`, `"vasp_std_2D"`, `"vasp_gam"`.
- **`qe`**: Name of the Quantum Espresso module, such as `"pw.x"`, `"ph.x"`.
- **`pwdft`**: Path to the PWDFT module.
- **`moduleload`**: Specifies the command for loading necessary modules, e.g., `load "vasp/5.4.4-intel2019,intelmpi/2019.update3,mkl/2019.update3"`. Separate modules with commas.
- **`genmode`**: Sets the mode for generating input files:
  - `"basic"`: Standard input file generation
  - `"vaspkit"`: Uses Vaspkit for generating input files
  - `"pmg"`: Uses pymatgen for input file generation
  - `"qe"`: Generates input files for Quantum Espresso
  - `"pwdft"`: Generates input files for PWDFT calculations
- **`kwargs`**: Additional parameters specific to the `genmode`.

**Use Cases**

```python
from pyhtstack2d.tools.genInput import GenRunDir
# Example: Generate input files in basic mode with custom INCAR settings and single-task submission script `run.sh`
user_incar_settings = {"NPAR": 16, "EDIFF": 1e-8, "IVDW": 12}
GenRunDir(posdir="TMD_POSCAR_dir",genmode="basic",inputOpt='scf',user_incar_settings=user_incar_settings).genInputfile()

# Example: Generate input files for multi-task calculations with Vaspkit and SLURM submission
GenRunDir(genmode="vaspkit",subset="slurm",kmeshrv=0.02,shoverwrite=True,taskname="scf").genInputfile()

# Example: Generate input files for a nested multi-level directory structure and multiple tasks ('relax', 'scf', and 'band')
GenRunDir(posdir="BiPOSCAR_dir",multilevel=4,inputOpt="relax-dip",kmeshrv=0.04,taskname="opt").genInputfile()
GenRunDir(posdir="BiPOSCAR_dir",multilevel=4,inputOpt="scf-dip",kmeshrv=0.04,taskname="scf").genInputfile()
GenRunDir(posdir="BiPOSCAR_dir",multilevel=4,inputOpt="band-dip",kmeshrv=0.04,postype='H',taskname="band").genInputfile()
```

### 4.3.7. `pyhtstack2d.calcSets.vaspMagSetsWriter`

This module generates the INCAR file for VASP calculations in magnetic systems.

**Parameters**

- **`maglist`**: Path to a file (TXT or JSON) specifying the location of magnetic systems. JSON files (e.g., `magdict.json`) are loaded as dictionaries for initializing `magdict`. Use `pyhtstack2d/Script/magScript` to generate `maglist.txt`.
- **`magpath`**: Path to the VASP calculation results for obtaining magnetic moments. For example, `magpath="scf/"` refers to the `OSZICAR` file in the `scf` directory. If `None`, magnetic moments are taken from `OSZICAR` in the main directory.
- **`magnetfile`**: Path to the file with magnetic moments. Ensure `LORBIT = 11` is set in the INCAR file.
- **`pmg`**: Enables the pymatgen module to retrieve magnetic moments.
- **`maglowerl`**: Lower limit for identifying magnetic atoms by magnetic moment.
- **`nomagatom`**: List of atoms that are typically non-magnetic (e.g., ["Cl", "F", "I", "Br"]), even if VASP magnetic projections show significant values on these atoms.
- **`natomlayer1`**: Number of atoms in the first layer.
- **`skipmagnetfile`**: If `True`, skips magnetic atom identification and directly sets specified atoms as magnetic.
- **`magatomlist`**: List of magnetic atoms provided directly when `skipmagnetfile` is `True`.
- **`magmom_setting`**: Sets specific magnetic moments for atoms, in the format `{"Fe": 3.0}`.
- **`pU`**: Adds a U correction to the magnetic atoms.
- **`onlymagpu`**: If `True`, applies U correction only to magnetic atoms.
- **`mkfm`**: Forces the creation of FM (ferromagnetic) folders.

For generating the INCAR file in magnetic system VASP calculations, use the following:

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspMagSetsWriter import MonoMagState
# Generate input files for magnetic monolayer calculation
# `INCARbasic` is the basic INCAR file. `supercellgen` is used to create supercells ("pmg" or "vaspkit").
# `vaspkitKpoints`, `kmseshScheme`, and `kmeshrv` are parameters for generating the KPOINTS file with Vaspkit.
MonoMagState().genInputfile(INCARbasic="INCAR-basic", supercellgen="pmg", vaspkitKpoints=True, kmseshScheme=2, kmeshrv=0.04)

from pyhtstack2d.calcSets.vaspMagSetsWriter import BiMagState
# Generate input files for bilayer magnetic calculations
# Please ensure that the `monomagdict.json` file is present and contains the magnetic information for the monolayer.
BiMagState(maglist="bimaterials.txt", pmg=True).genInputfile(INCARbasic="INCAR")
```

> **Note**: The file `biInputdict.json` is automatically generated when `BiMagState()` is called. If the magnetic atom or monolayer information changes, delete `biInputdict.json` before calling `BiMagState()` again.


The keys in `biInputdict.json` are described as follows:

| Key             | Description                                                                                                    |
|-----------------|----------------------------------------------------------------------------------------------------------------|
| `Inputpath`     | Primary path to input file settings                                                                           |
| `magatom`       | List of magnetic atoms                                                                                        |
| `magmomlist`    | List of magnetic moments                                                                                      |
| `submaterial`   | Subpaths that share input files with the primary path                                                         |
| `natomlayer1`   | Number of atoms in the first layer                                                                            |
| `blabel`        | Label of the bilayer structure                                                                                |
| `bl1label`      | Label of the first layer                                                                                      |
| `bl2label`      | Label of the second layer                                                                                     |
| `element`       | Element list from the POSCAR file                                                                             |
| `nmag`          | Number of magnetic atoms                                                                                      |
| `magstate`      | Magnetic state of the bilayer structure (e.g., `33` for FM/FM). Codes: `0` - NM, `1` - AFM1, `2` - AFM2, `3` - FM, `4` - AFM3, `5` - unknown |
| `atommagstate`  | Magnetic state of the layer where each atom is located                                                        |
| `atomlayer`     | The specific layer where each atom is located                                                                 |



Since VASP performs ferromagnetic (FM) calculations by default when spin-polarization is enabled, FM calculations can be optionally skipped in practical cases to reduce computational costs.

Job submission scripts can be generated with the function `genRunsh`.

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspMagSetsWriter import MonoMagState
# Generate the PBS job submission script for magnetic monolayer calculations
MonoMagState().genRunsh(pbs=True, moduleload="vasp/5.4.4-intel2019", vasp="vasp_std_2D")

from pyhtstack2d.calcSets.vaspMagSetsWriter import BiMagState
# Generate the SLURM job submission script for bilayer magnetic calculations
BiMagState().genRunsh(pbs=True, moduleload="vasp/5.4.4-intel2019", vasp="vasp_std_2D")
```

To obtain results such as analyzing the magnetic ground state (MGS), use the following code:

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspMagSetsWriter import MonoMagState
# Parameters:
# - oszipath: Path to the FM calculation results.
# - FMm4: If True, multiplies the FM calculation energy by 4.
# - copyCONTCAR: If True, copies the CONTCAR file of the MGS to the POSCAR file.
MonoMagState(maglist="magdict.json").getEnergy(oszipath=None, FMm4=False, copyCONTCAR=True)

from pyhtstack2d.calcSets.vaspMagSetsWriter import BiMagState
# - skipNFM: If True, skips NM/FM calculations for bilayer structures. However, it will still read energies for NM/NM and NM/FM combinations if set to True.
BiMagState().getEnergy(oszipath=None, FMm4=False, copyCONTCAR=True, skipNFM=False)
```

The results are saved as JSON files:
- `MonoMagState().getEnergy()` results are saved to `"monomagdict.json"`.
- `BiMagState().getEnergy()` results are saved to `"bimagdict.json"` and `"bimagdictSum.json"`.

The keys and associated descriptions of `monomagdict.json` are as follows:

| Key         | Description                                                                            |
|-------------|----------------------------------------------------------------------------------------|
| `material`  | Path to the calculation for the material                                               |
| `magatom`   | Magnetic atom                                                                          |
| `magE`      | Energy of the magnetic calculation, covering FM and Antiferromagnetic (AFM) states     |
| `magstate`  | The magnetic ground state (MGS)                                                        |

The keys and associated descriptions of `bimagdict.json` or `bimagdictSum.json` are as follows:

| Key                  | Description                                                                                              |
|----------------------|----------------------------------------------------------------------------------------------------------|
| `stackmode`          | Stacking mode of the bilayer structure                                                                   |
| `magstatelist`       | MGS of the bilayer structure for each stacking mode                                                      |
| `magstatelistcorrect`| Corrected MGS after re-identifying magnetic atoms and their moments                                      |
| `maglist`            | Total magnetic moment of the MGS                                                                         |
| `magE`               | Energy of the magnetic calculation, including NM/FM and AFM states                                       |
| `mag`                | Total magnetic moment of the magnetic calculation, covering NM/FM and AFM states                         |
| `NFMstate`           | State of the NM/FM calculation                                                                           |
| `mixstatecorrect`    | Indicates if the magnetic state is mixed with the stacking mode                                          |


### 4.3.8. `pyhtstack2d.calcSets.Pu`

This module adds U corrections to the INCAR file based on a given POSCAR file.

**Parameters**

- **`elemlist`**: List of elements in the structure.
- **`incar`**: Path to the INCAR file. If provided, U settings will be directly added to this file. Otherwise, the function returns a dictionary with the U settings.
- **`u_setting`**: Dictionary specifying U values for particular elements.
- **`magatomlist`**: List of magnetic atoms. If specified, U corrections will only be applied to these atoms.
- **`openmixparam`**: If `True`, enables the mixing parameters in the INCAR file.
- **`mixparam`**: Specifies the mixing parameters.

**Use Case**

```python
from pyhtstack2d.calcSets.Pu import INCARPu

# Example: Adding U corrections to elements in the INCAR file
elem = "Cl   Bi   S    Se   Sc".split()
INCARPu(elem, "INCAR")
```


### 4.3.9. `pyhtstack2d.calcSets.vaspINCARupdate`

This module updates INCAR files with new parameters or settings from an `INCAR-basic` file.

**Parameters**

- **`materiallist`**: List of materials.
- **`dirprefix`**: Prefix for directories containing materials (e.g., `"scf"`).
- **`incarbasic`**: Path to the basic INCAR file or a dictionary of parameters.
- **`pU`**: If `True`, adds U settings to the INCAR file.
- **`collateincar`**: If `True`, collates INCAR content after updating to remove duplicate lines and rearrange parameters.
- **`ismag`**: Specifies if the calculation is magnetic. If `True`, it will search for directories such as `FM`, `AFM1`, `AFM2`, etc.

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspINCARupdate import UpdateINCAR
# Example: Update all INCAR files in specified directories using parameters from `incarbasic="INCAR"` file.
UpdateINCAR(materiallist="maglist.txt", incarbasic="INCAR" ,pU=True, ismag=True).updateIncarfiles()


# Example: Directly overwrite the original INCAR file and collate content to remove duplicates.
UpdateINCAR(materiallist="maglist.txt", collateincar=True).genIncarfiles()
```


## 4.4. Generating Job Submission Scripts

### 4.4.1 `pyhtstack2d.tools.genInput.GenMultiSh`

This module generates job submission scripts for SLURM or PBS job schedulers.

In **4.3.6** (`pyhtstack2d.tools.genInput.GenRunDir`), calling the `genInputfile()` function with `taskname=None` generates input files and a single-task job submission script (`run.sh`). For multi-task job submissions, use the `GenMultiSh` class in the `pyhtstack2d.tools.genInput` module.

**Parameters**

- **`tasklist`**: The list of task names, e.g., `["relax", "scf", "band"]`.
- **`workdir`**: The directory where the `run.sh` file is generated.
- **`multilevel`**: Same as the `multilevel` parameter in the `GenRunDir()` class; specifies the level of nested directories.
- **`incpath`**, **`kppath`**, **`potpath`**: Paths to shared INCAR, KPOINTS, or POTCAR files for all calculations if these parameters are not `None`.
- **`opts`**: If the task involves structural optimization, a secondary structural optimization calculation is automatically triggered to avoid local minima in the first optimization.
- **`failact`**: Specifies the action on calculation failure. `"continue"`: Continues to the next task calculation; `"break"`: Stops remaining tasks and moves to other material calculations.
- **`saveall`**: If `True`, saves all output files for each task in corresponding directories. If `False`, only the files specified in `savefiles` are saved.
- **`savefiles`**: The list of output files to save, e.g., `["OUTCAR", "CONTCAR"]`.
- **`checkf`**: Checks that the input files are completely prepared before execution.

The parameters **`multilevel`**, **`subset`**, **`mpirun`**, **`vasp`**, and **`moduleload`** function the same as in the `GenRunDir()` class.

**Use Cases**

```python
from pyhtstack2d.tools.genInput import GenMultiSh
# Example: Generate job submission script with all output files saved for tasks "relax", "scf", and "band"
GenMultiSh(tasklist=["relax", "scf", "band"], saveall=True).genRunsh()

# Example: When INCAR and KPOINTS files are shared across all materials
incpath = ["Inputfiles/INCAR-relax", "Inputfiles/INCAR-scf", "Inputfiles/INCAR-band"]
kppath = ["Inputfiles/KPOINTS-relax", "Inputfiles/KPOINTS-scf", "Inputfiles/KPOINTS-band"]
GenMultiSh(tasklist=["relax", "scf", "band"], incpath=incpath, kppath=kppath).genRunsh()

# Example: For Quantum Espresso (QE) multi-task calculations
from pyhtstack2d.tools.genInput import GenMultiQESh
GenMultiQESh(tasklist=["relax", "scf", "band"]).genRunsh()
```

When encountering structural calculation failures, such as non-convergence, you can modify parameter settings and generate job submission scripts to rerun the calculations for specific paths. The following code demonstrates how to generate job submission scripts for rerunning the structures under specified paths:

```python
from pyhtstack2d.tools.genInput import GenRunDir
# For single-task job submission rerun
GenRunDir().failrw()

# For multi-task job submission in nested directories
from pyhtstack2d.tools.genInput import MulTaskRW
MulTaskRW(multilevel=4, saveall=True).genRunsh()
```

## 4.5. Batch Extraction and Simple Analysis of Results

### 4.5.1. `pyhtstack2d.tools.genInput.OptPOSCAR`

This module extracts optimized structures from `CONTCAR` files by generating a `getPOSCAR.sh` script.

**Parameters**

- **`workdir`**: Path to the working directory. If not specified, the current directory (`os.getcwd()`) is used.
- **`posdir`**: Directory to save the optimized `POSCAR` files.
- **`optpath`**: Path to the directory containing optimization calculation results. For multi-task workflows, `CONTCAR` files may be located in an `opt/` subdirectory.
- **`multilevel`**: Specifies the level of nested directories to consider.

**Use Case**

```python
from pyhtstack2d.tools.genInput import OptPOSCAR
# Generate the getPOSCAR.sh script for batch extraction of optimized structures
OptPOSCAR().genGetsh()
```

### 4.5.2. `pyhtstack2d.tools.genInput.GetMagEntropy`

This module retrieves magnetic moments and entropy values after VASP calculations. It is particularly useful for separating materials with mixed magnetic and non-magnetic properties. If a material's entropy exceeds a certain threshold, it may require recalculations to ensure that the calculated entropy per atom (`T*S`) is below 1 meV.

**Parameters**

- **`workdir`**: Path to the working directory. Defaults to the current directory (`os.getcwd()`).
- **`magpath`**: Path to the VASP calculation results for retrieving magnetic moments. For example, `magpath="scf/"` indicates that magnetic moments are in the `OSZICAR` file within the `scf` directory; `magpath="-scf"` refers to `OSZICAR-scf` in the material directory.
- **`multilevel`**: Specifies the level of nested directories to consider.

**Use Cases**

```python
from pyhtstack2d.tools.genInput import GetMagEntropy
# Example: Separate materials based on entropy values
# - `entxt`: Text file containing entropy data. If not present, entropy information will be extracted from `workdir`.
# - `entropytol`: Threshold for entropy; materials with higher entropy will be marked for recalculation.
# - `recalheader`: The path of header for generating job submission scripts, used when recalculating materials.
GetMagEntropy(magpath="-opt").separateEntropy(entxt="entropy.txt", entropytol=0.001, recalheader="header.sh")

# Example: Separate materials based on magnetic moments
# - `magtxt`: Text file with total magnetic moment data. If not present, magnetic moment data will be extracted from `workdir`.
# - `magtol`: Magnetic moment threshold; materials exceeding this threshold will be added to the magnetic materials list.
# This function generates `maglist.txt` and `nonmaglist.txt` in `workdir`.
GetMagEntropy(magpath="opt/").separateMag(magtxt="mag.txt", magtol=0.5)
```


### 4.5.3. `pyhtstack2d.analysis.extractResults`

This module is designed for batch extraction and analysis of VASP calculation results.

**Parameters**

- **`materialfile`**: Path to a file listing materials for extraction. If not provided, the module will traverse the current folder.
- **`workdir`**: Path to the working directory. Defaults to the current directory.
- **`multilevel`**: Specifies the depth of nested directories to consider.
- **`scf`**: Name of the file or folder containing the SCF or optimization/relaxation results.
- **`band`**, **`hybband`**, **`optical`**, **`eps`**: Names of files or folders containing results for these respective calculations.
- **`directbandgap`**: Extracts direct band gap information for indirect gap materials.
- **`bandplot`**: Determines whether to plot the band structure or layer-projected band structure.
- **`banderange`**: Specifies the energy range for band structure plotting.
- **`imag_format`**: Sets the image format for band structure plots, e.g., `"png"` or `"pdf"`.
- **`elempro`**: Extracts element-projected band structure/DOS information.
- **`dim`**: Dimension of the system (1D, 2D, or 3D).
- **`energyunit`**: Sets the energy unit for the absorption coefficient. Options are:
  - `1`: eV
  - `2`: nm
  - `3`: THz
  - `4`: cm⁻¹
- **`weightsave`**: Saves the weights extracted from the PROCAR file if `True`.
- **`nlayer`**: Number of layers for 2D materials.
- **`mag`**: Extracts magnetic moment information for the system.
- **`pmg`**: Uses pymatgen to extract calculation results if `True`.
- **`vaspkit`**: Uses Vaspkit to extract calculation results if `True`.
- **`infodict`**: Loads JSON content directly for analysis, bypassing batch extraction.

> When calling this class without the `infodict` parameter, an `"info.json"` file is automatically generated. If `"info.json"` already exists, the file will be saved as `"info1.json"`, `"info2.json"`, and so forth.

**Use Cases**

```python
from pyhtstack2d.analysis.extractResults import GetResults
# For monolayers.
GetResults(scf="scf", band="band", pmg=False, multilevel=4, mag=True)

# For bilayers:
# Calling `getInfoBi()` will analyze the binding energy and band alignment type.
# The `getInfoBi()` function requires `monoinfo.json`, which contains monolayer information, to be in the current directory.
# This function generates `biinfo.json` and `biinfo-simple.json`.
GetResults(multilevel=4, nlayer=2, infodict="info.json").getInfoBi(monodict="monoinfo.json")
```

The keys and associated descriptions of `info.json` are as follows:

| Key                  | Description                                                                                        |
|----------------------|----------------------------------------------------------------------------------------------------|
| `elements`           | List of elements                                                                                   |
| `label`              | Label for the material                                                                             |
| `atom_indices`       | Atomic height ordering in the z-direction                                                          |
| `lattice`            | Lattice matrix                                                                                     |
| `abc`                | Lattice constants                                                                                  |
| `cart_coords`        | Atomic positions in Cartesian coordinates                                                          |
| `zdiff`              | Atomic spacing along the z-axis                                                                    |
| `spacegroup`         | Space group of the structure                                                                       |
| `Efermi`             | Fermi level energy                                                                                |
| `E0`                 | Total energy                                                                                       |
| `mag`                | Magnetic moment                                                                                    |
| `Magnetization_x`    | Atomic projected magnetic moment                                                                   |
| `bandgap`            | Minimum bandgap                                                                                   |
| `direct`             | Indicates if the bandgap is direct (`True`) or indirect (`False`)                                 |
| `is_metals`          | Indicates if the system is metallic (`True` or `False`)                                           |
| `vbm` or `cbm`       | Energy of the valence band maximum (VBM) or conduction band minimum (CBM)                         |
| `vbm_index` or `cbm_index` | Index information for VBM or CBM, including spin direction, band index, and k-point index  |
| `vbmpro` or `cbmpro` | Atomic projection contribution at the VBM or CBM                                                   |
| `directbandgap`      | Contains direct bandgap information for spins 1 and -1, available when `pmg` is `True`             |
| `vacuum_potential`   | Vacuum potential energy                                                                            |
| `work_function`      | Work function energy                                                                               |
| `coords_z`           | Grid along the z-direction for plotting                                                            |
| `potential_z`        | Potential energies along the z-direction                                                           |
| `vacuum_potential_dict` | Dictionary containing left, right, and differential vacuum potential energy                   |
| `absorption`         | Absorption data                                                                                    |


The special keys and associated descriptions of `biinfo.json` are as follows:

| Key               | Description                                                                                          |
|-------------------|------------------------------------------------------------------------------------------------------|
| `mismatch`        | Lattice mismatch between layers                                                                      |
| `natomlayer1`     | Number of atoms in the first layer                                                                   |
| `blabel`          | Label for the bilayer                                                                                |
| `bl1label` or `bl2label` | Label for the first or second layer                                                          |
| `l1info` or `l2info` | Information for the first or second layer                                                        |
| `stackmode`       | Stacking mode of the bilayer                                                                         |
| `Eb`              | Binding energy of the bilayer                                                                        |
| `dinterlayer`     | Interlayer spacing between the two layers                                                            |
| `vbmprosum` or `cbmprosum` | Layer-summed atomic projection contributions at the VBM or CBM                              |
| `bandalignment`   | Type of band alignment                                                                               |
| `deltabandgap`    | Difference between the minimum bandgap and direct bandgap                                            |

In `biinfo-simple.json`, single-layer information and elemental details are omitted, providing a concise display of bandgaps, binding energies, and band alignment types for different stacking configurations.


### 4.5.4. `pyhtstack2d.analysis.plotBAND`

This module is used to plot the band structure.

**Parameters**

- **`bandpath`**: Path to the band structure file.
- **`indices`**: Indices of the atoms in the layer.
- **`natomlayer1`**: Number of atoms in the first layer.
- **`hybrid`**: Whether to plot the hybrid band structure.
- **`erange`**: Energy range for band structure plotting.
- **`elempro`**: Whether to plot the element-projected band structure.
- **`pmg`**: If `True`, uses pymatgen to plot the band structure.
- **`BSDOSPlotter`**: `BSDOSPlotter` module from pymatgen for plotting.
- **`BSPlotter`**: `BSPlotter` module from pymatgen for plotting.
- **`Vasprun`**: `Vasprun` module from pymatgen to process VASP results.
- **`vaspkit`**: If `True`, uses Vaspkit to generate the `band.dat` file.
- **`imag_format`**: Image format for the saved plot (e.g., `"png"`, `"pdf"`).
- **`dpi`**: Resolution (dots per inch) for the saved plot.

**Use Cases**

```python
from pyhtstack2d.analysis.plotBAND import plotBS
import os
# Example: Using Vaspkit to extract `band.dat` and plot band structure
plotBS(bandpath=os.getcwd(), pmg=False).plotbsdos()

# Example: Using pymatgen's Vasprun and BSDOSPlotter modules to plot the band structure
Vasprun = __import__('pymatgen.io.vasp.outputs', fromlist=['Vasprun']).Vasprun
BSDOSPlotter = __import__('pymatgen.electronic_structure.plotter', fromlist=['BSDOSPlotter']).BSDOSPlotter
plotBS(bandpath=os.getcwd(), pmg=True, Vasprun=Vasprun, BSDOSPlotter=BSDOSPlotter).plotbsdos()

# Example: Plotting a layer-projected band structure for specified atom indices
plotBS(bandpath=os.getcwd(), pmg=False, natomlayer1=3, indices=[0, 1, 2, 3, 4, 5]).plotbsdos()
```


# 5. High-throughput illustrative examples

## 5.1. Hexagonal crystal systems containing 3 atoms monolayers

For more illustrative examples, see Examples 1-3 covered in our published papers.

Below, we provide an additional example, demonstrating how to identify monolayers by their band types, and subsequently perform high-throughput stacking to identify combinations that facilitate an indirect-to-direct bandgap transition.

```python
# Import necessary modules
from pyhtstack2d.tools.queryDB import getuid
import pandas

# Define database path and selection criteria for materials
database = "../c2db.db"
criteria_p3m1 = 'natoms=3, dyn_stab=Yes, hform<0, gap>0, layergroup=p3m1, magstate=NM'
getuid(database, criteria_p3m1, save_csv_path='tmd_p3m1.csv', save_csv=True)

# Load data and separate materials based on direct and indirect bandgap
tmd_p3m1 = pandas.read_csv('tmd_p3m1.csv')
tmd_p3m1_uid_indir = []
for i in range(len(tmd_p3m1)):
    if abs(tmd_p3m1.iloc[i]['gap'] - tmd_p3m1.iloc[i]['gap_dir']) > 1e-3:
        tmd_p3m1_uid_indir.append(tmd_p3m1.iloc[i]['uid'])

# Fetch POSCAR files for selected UIDs
from pyhtstack2d.tools.queryDB import getPOSACRfromuid
getPOSACRfromuid(database, tmd_p3m1_uid_indir, overwrite=False, save_path='POSCAR_Gind')

# Generate bilayers from indirect bandgap materials
from pyhtstack2d.buildbilayer.batchStackBilayer import GenBiLayer
GenBiLayer('POSCAR_Gind', genmode="tmdh", la_mismatch=5.0, homo=True).batch_stack()

# Use vaspkit to prepare POTCAR and KPOINTS for each POSCAR.
from pyhtstack2d.tools.genInput import GenRunDir, GenMultiSh
GenRunDir(genmode="vaspkit", posdir="BiPOSCAR_dir", workdir="BiPOSCAR_dir", multilevel=4, taskname="opt", incexis=True).genInputfile()

# Generate batch job submission scripts for optimization, self-consistent field, and band structure calculations
# Define paths to the INCAR files shared across different types of calculations
incpath = ["INPUT_file/INCAR-opt", "INPUT_file/INCAR-scf", "INPUT_file/INCAR-band"]
GenMultiSh(tasklist=["opt", "scf", "band"], workdir="BiPOSCAR_dir", multilevel=4, opts=True, saveall=True, incpath=incpath).genRunsh()
```

## 5.2. Homocrystalline system with multiple atomic number monolayers

The above examples all limit the need for all monolayers to have the same number of atoms and crystal system, and here an example is given that allows stacking between monolayers of the same crystal system, but without limiting the number of atoms.

To be continued and will be added soon.

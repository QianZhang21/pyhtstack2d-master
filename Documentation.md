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

```

### 4.1.2 `pyhtstack2d.tools.queryDBurl`

This module retrieves the POSCAR file for a specified material directly from a webpage URL using its unique identifier (uid).

**Use cases**:

```python
from pyhtstack2d.tools.getPOSCARurl import DBuidPOSCAR

```

### 4.1.3 `pyhtstack2d.tools.cif2pos`

This module facilitates the conversion of CIF files to POSCAR format.

**Use cases**:

```python

```

### 4.1.3 `pyhtstack2d.tools.POSCARelemc`

This module is used to modify the elements in a POSCAR file.

**Use cases**:

```python

```

### 4.1.4 `pyhtstack2d.tools.renameposcar`

For stacking purposes, POSCAR filenames need to follow the "uid-POSCAR" format. This module renames POSCAR files that are in the "POSCAR-uid" format.

**Use cases**:

```python

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

- **`savenamemode`**: Naming mode for saved POSCAR files.
  - **1**: Files saved in `"{savepath}/{formula}_{mismatch}/{genmode}/{cord*}/POSCAR"`.
  - **2**: Files saved as `"{savepath}/{formula}_{genmode}_{cord*}-{mismatch}-POSCAR"`.
  
**Use Cases:**

For instance, you can stack MoS₂ and WSe₂ bilayer structures by specifying the elements, lattice parameters, and atomic positions, or by using POSCAR files.

```python
from pyhtstack2d.buildbilayer.stackBilayer import Bilayer

```

### 4.2.2. `pyhtstack2d.buildbilayer.batchStackBilayer`

This module provides a class for batch stacking of bilayer structures.

**Parameters:**

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

```

## 4.3. Generating Input Files

### Fundamental Parameters

- **`taskname`**: Name used as a suffix or prefix for input filenames in multi-task calculations.

These fundamental parameters are referenced across modules 4.3.1 to 4.3.5.

### 4.3.1. `pyhtstack2d.calcSets.qeSets`

This module generates input files for Quantum Espresso (QE) calculations based on POSCAR files.

**Special Parameters**

- **`spinOpen`**: Indicates if the input file should be configured for spin-polarized calculations.

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

- **`is_print`**: If enabled, prints the path where the input files are saved.

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspSetsWriter import InputWriter
import os
from glob import glob

```

### 4.3.4. `pyhtstack2d.calcSets.vaspPmgWriter`

This module generates input files for VASP calculations using the pymatgen interface, based on POSCAR files.

**Special Parameters**

- **`is_print`**: If enabled, prints the path where the input files are saved.

**Use Cases**

```python
from pyhtstack2d.calcSets.vaspPmgWriter import PmgInputWriter

```

### 4.3.5. `pyhtstack2d.calcSets.vaspVaspkitWriter`

This module generates input files for VASP calculations using the Vaspkit interface, based on POSCAR files. Typically, this module creates a `mkdir.sh` script that needs to be executed to complete the full setup of input files.

**Special Parameters**

- **`is_print`**: If enabled, prints the path where input files are saved.

**Use Cases**

```python

```

### 4.3.6. `pyhtstack2d.tools.genInput.GenRunDir`

This module is designed for batch generation of input files for VASP calculations, with flexible settings for handling various submission and directory structures.

**Parameters**

- **`kwargs`**: Additional parameters specific to the `genmode`.

**Use Cases**

```python

```

### 4.3.7. `pyhtstack2d.calcSets.vaspMagSetsWriter`

This module generates the INCAR file for VASP calculations in magnetic systems.

**Parameters**

- **`mkfm`**: Forces the creation of FM (ferromagnetic) folders.

For generating the INCAR file in magnetic system VASP calculations, use the following:

**Use Cases**

```python

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

```

To obtain results such as analyzing the magnetic ground state (MGS), use the following code:

**Use Cases**

```python

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

- **`mixparam`**: Specifies the mixing parameters.

**Use Case**

```python

```


### 4.3.9. `pyhtstack2d.calcSets.vaspINCARupdate`

This module updates INCAR files with new parameters or settings from an `INCAR-basic` file.

**Parameters**

- **`ismag`**: Specifies if the calculation is magnetic. If `True`, it will search for directories such as `FM`, `AFM1`, `AFM2`, etc.

**Use Cases**

```python

```


## 4.4. Generating Job Submission Scripts

### 4.4.1 `pyhtstack2d.tools.genInput.GenMultiSh`

This module generates job submission scripts for SLURM or PBS job schedulers.

In **4.3.6** (`pyhtstack2d.tools.genInput.GenRunDir`), calling the `genInputfile()` function with `taskname=None` generates input files and a single-task job submission script (`run.sh`). For multi-task job submissions, use the `GenMultiSh` class in the `pyhtstack2d.tools.genInput` module.

**Parameters**

- **`checkf`**: Checks that the input files are completely prepared before execution.

The parameters **`multilevel`**, **`subset`**, **`mpirun`**, **`vasp`**, and **`moduleload`** function the same as in the `GenRunDir()` class.

**Use Cases**

```python

```

When encountering structural calculation failures, such as non-convergence, you can modify parameter settings and generate job submission scripts to rerun the calculations for specific paths. The following code demonstrates how to generate job submission scripts for rerunning the structures under specified paths:

```python

```

## 4.5. Batch Extraction and Simple Analysis of Results

### 4.5.1. `pyhtstack2d.tools.genInput.OptPOSCAR`

This module extracts optimized structures from `CONTCAR` files by generating a `getPOSCAR.sh` script.

**Parameters**

- **`multilevel`**: Specifies the level of nested directories to consider.

**Use Case**

```python

```

### 4.5.2. `pyhtstack2d.tools.genInput.GetMagEntropy`

This module retrieves magnetic moments and entropy values after VASP calculations. It is particularly useful for separating materials with mixed magnetic and non-magnetic properties. If a material's entropy exceeds a certain threshold, it may require recalculations to ensure that the calculated entropy per atom (`T*S`) is below 1 meV.

**Parameters**

- **`multilevel`**: Specifies the level of nested directories to consider.

**Use Cases**

```python

```


### 4.5.3. `pyhtstack2d.analysis.extractResults`

This module is designed for batch extraction and analysis of VASP calculation results.

**Parameters**

- **`infodict`**: Loads JSON content directly for analysis, bypassing batch extraction.

> When calling this class without the `infodict` parameter, an `"info.json"` file is automatically generated. If `"info.json"` already exists, the file will be saved as `"info1.json"`, `"info2.json"`, and so forth.

**Use Cases**

```python

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

- **`dpi`**: Resolution (dots per inch) for the saved plot.

**Use Cases**

```python

```


# 5. High-throughput illustrative examples

For more illustrative examples, see Examples 1-3 covered in our published papers.

Below, we provide an additional example, demonstrating how to identify monolayers by their band types, and subsequently perform high-throughput stacking to identify combinations that facilitate an indirect-to-direct bandgap transition.

```python

```


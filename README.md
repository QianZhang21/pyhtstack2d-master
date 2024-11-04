# PyHTStack2D: a Python package for high-throughput homo/hetero stacking of 2D materials

## Brief Introduction
PyHTStack2D is a Python package designed for high-throughput stacking of two-dimensional (2D) layered materials. It assists in generating high-throughput job submission scripts and batch extraction of results.

## Features
- **POSCAR Retrieval:** Fetch POSCAR files of 2D monolayer materials from the C2DB database.
- **Layer Stacking:** Create bilayer structures through stacking, considering multiple high-symmetry stacking modes for special systems. When any layer exhibits mirror asymmetry, it considers flipping it to generate new bilayer structures with different stacking orders.
- **Input File Generation:** Based on the POSCAR files and specified parameters, generate the necessary input files for first-principles calculations.
- **Script Generation:** Create effective shell scripts for batch job submissions.
- **Data Extraction and Analysis:** Perform batch extraction of calculation results and conduct basic analyses.

## Installation and Use
For detailed installation instructions, configuration settings, and usage examples, please refer to the `Documentation.md` file included with the package or visit our [GitHub repository](https://github.com/QianZhang21/PyHTStack2D).

## License

This software is licensed under the BSD 3-Clause License, one of the more permissive free software licenses. This license allows you to use, modify, and distribute the software either in source code or binary form. For the specific terms and conditions, refer to the full license text available online at [opensource.org](https://opensource.org/licenses/BSD-3-Clause).

## How to Cite
To cite PyHTStack2D in your research, please use the following citation:

[1] Qian Zhang, Wei Hu, and Jinlong Yang, "PyHTStack2D: a Python package for high-throughput homo/hetero stacking of 2D materials," 


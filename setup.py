import os
import json
from setuptools import setup, find_packages
from setuptools.command.install import install


class CustomInstall(install):
    def run(self):
        # The default configuration values
        config_dict_set = {
            "POTCAR_PATH": "Specify your path here",
            "PBS_set": {
                "PBS -N ": "jobname",
                "PBS -l ": "nodes=1:ppn=28",
                "PBS -j ": "oe",
                "PBS -q ": "normal"
            },
            "PBS_echo": [
                "echo $PBS_O_WORKDIR",
                "NCPUS=`wc -l $PBS_NODEFILE | awk '{print $1}'`",
                "echo $NCPUS"
            ],
            "SLURM_set": {
                "SBATCH -J ": "jobname",
                "SBATCH -N ": "1",
                "SBATCH -n ": "48",
                "SBATCH -p ": "normal",
                "SBATCH -o ": "%j.log",
                "SBATCH -e ": "%j.err"
            },
            "SLURM_echo": [
                "echo Time is `date`",
                "echo Directory is $PWD",
                "echo This job runs on the following nodes:",
                "echo $SLURM_JOB_NODELIST",
                "echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores."
            ]
        }

        print("\nDefault configuration values:")
        for key, value in config_dict_set.items():
            if key in ["PBS_echo", "SLURM_echo"]:
                continue
            print(f"\n{key}: {value}")
            if isinstance(value, dict):
                for subkey, subvalue in value.items():
                    new_value = input(f"Enter value for '{subkey}' (default: {subvalue}): ") or subvalue
                    config_dict_set[key][subkey] = new_value
            else:
                new_value = input(f"Enter value for '{key}' (default: {value}): ") or value
                config_dict_set[key] = new_value

        # The path to the configuration file
        config_dir = os.path.join(os.path.expanduser("~"), ".config")
        config_path = os.path.join(config_dir, ".PyHTStack2D.json")

        # Ensure the directory exists
        os.makedirs(config_dir, exist_ok=True)

        # Save the configuration to the file
        with open(config_path, 'w') as f:
            json.dump(config_dict_set, f, indent=4)

        print(f"\nConfiguration saved at {config_path}")

        # Ask the user if they want to install the recommended dependencies
        choice = input("Do you want to install recommended dependencies? (yes/no): ").strip().lower()
        if choice in ['yes', 'y']:
            # Install the recommended dependencies
            os.system("pip install ase==3.22.1 matplotlib==3.5.1 numpy==2.1.2 pandas==1.4.4 pymatgen==2023.11.12")
            print("\nRecommended dependencies installed.")
        else:
            print("\nSkipping recommended dependencies installation.")
        # Continue with the installation
        install.run(self)


# Configure the package
setup(
    name="pyhtstack2d",
    version="0.1",
    author="Qian Zhang, Wei Hu, Jinlong Yang, University of Science and Technology of China, Hefei, China.",
    license="BSD-3-Clause",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "matplotlib",
        "ase",
        "pandas",
        "pymatgen"
    ],
    extras_require={
        "recommended": [  # Recommended dependencies version
            "ase==3.24.0",
            "matplotlib==3.10.1",
            "numpy==2.2.4",
            "pandas==2.2.3",
            "pymatgen==2024.10.3"
        ]},
    cmdclass={'install': CustomInstall},  # Use the custom install class
)

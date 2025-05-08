import os
import json
import subprocess


def configure_and_install():
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

        if isinstance(value, dict):
            for subkey, subvalue in value.items():
                new_value = input(f"Enter value for '{subkey}' (default: {subvalue}): ") or subvalue
                config_dict_set[key][subkey] = new_value
        else:
            new_value = input(f"Enter value for '{key}' (default: {value}): ") or value
            config_dict_set[key] = new_value

    config_dir = os.path.join(os.path.expanduser("~"), ".config")
    config_path = os.path.join(config_dir, ".PyHTStack2D.json")
    os.makedirs(config_dir, exist_ok=True)
    with open(config_path, 'w') as f:
        json.dump(config_dict_set, f, indent=4)
    print(f"\nConfiguration saved at {config_path}")

    choice = input("Do you want to install recommended dependencies? (yes/no): ").strip().lower()
    if choice in ['yes', 'y']:
        subprocess.run(["pip", "install", "ase==3.24.0", "matplotlib==3.10.1", "numpy==2.2.4", "pandas==2.2.3", "pymatgen==2024.10.3"])
        print("\nRecommended dependencies installed.")
    else:
        print("\nSkipping recommended dependencies installation.")

    subprocess.run(["pip", "install", "./dist/pyhtstack2d-0.1-py3-none-any.whl"])  # The location of the WHL file is subject to change on a case-by-case basis.


if __name__ == "__main__":
    configure_and_install()

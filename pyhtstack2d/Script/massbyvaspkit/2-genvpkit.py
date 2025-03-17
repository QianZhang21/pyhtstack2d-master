# Generate VPKIT.in for hexagonal 2D non-magnetic materials.
import os
import numpy as np
import subprocess


with open("mass_mid.txt", 'r') as f:
    midlist = f.readlines()


def run_subprocess(command):
    return subprocess.run(command, shell=True, capture_output=True, text=True).stdout.strip()


def get_band_location(bandpath, band_type, index, suffix=' (TO)'):
    command = f"grep 'Location of {band_type}{suffix}' {bandpath}/BAND_GAP | awk '{{print ${index}}}'"
    return float(run_subprocess(command))


def get_vbm_cbm_locations(bandpath, suffix=' (TO)'):
    if suffix:
        rangestar = 5
    else:
        rangestar = 4
    vbmlocation = [get_band_location(bandpath, "VBM", i, suffix) for i in range(rangestar, rangestar + 3)]
    cbmlocation = [get_band_location(bandpath, "CBM", i, suffix) for i in range(rangestar, rangestar + 3)]
    return vbmlocation, cbmlocation


Gpath = [0.000, 0.000, 0.000]
Kpath = [0.3333333333, 0.3333333333, 0.000]
Mpath = [0.500, 0.000, 0.000]


def writevpkit_kp(vpkitpath, location, label="V"):
    if np.linalg.norm([location[i] - Gpath[i] for i in range(3)]) < 0.01:
        with open(vpkitpath, "a", newline="\n") as f:
            f.write("0.000 0.000 0.000 0.500 0.000 0.000 Γ->M\n")
            f.write("0.000 0.000 0.000 0.3333333333 0.3333333333 0.000 Γ->K\n")
    elif np.linalg.norm([location[i] - Kpath[i] for i in range(3)]) < 0.01:
        with open(vpkitpath, "a", newline="\n") as f:
            f.write("0.3333333333 0.3333333333 0.000 0.500 0.000 0.000 K->M\n")
            f.write("0.3333333333 0.3333333333 0.000 0.000 0.000 0.000 K->Γ\n")
    elif np.linalg.norm([location[i] - Mpath[i] for i in range(3)]) < 0.01:
        with open(vpkitpath, "a", newline="\n") as f:
            f.write("0.500 0.000 0.000 0.3333333333 0.3333333333 0.000 M->K\n")
            f.write("0.500 0.000 0.000 0.000 0.000 0.000 M->Γ\n")
    else:
        if abs(location[1]-0.0) < 0.01 and 0.0 <= location[0] <= 0.5:
            with open(vpkitpath, "a", newline="\n") as f:
                f.write(f"{location[0]} {location[1]} {location[2]} 0.000 0.000 0.000 {label}->Γ\n")
                f.write(f"{location[0]} {location[1]} {location[2]} 0.500 0.000 0.000 {label}->M\n")
        elif 0.0 <= location[1] <= 1/3:
            if location[0] <= 1/3:
                with open(vpkitpath, "a", newline="\n") as f:
                    f.write(f"{location[0]} {location[1]} {location[2]} 0.3333333333 0.3333333333 0.000 {label}->K\n")
                    f.write(f"{location[0]} {location[1]} {location[2]} 0.000 0.000 0.000 {label}->Γ\n")
            else:
                with open(vpkitpath, "a", newline="\n") as f:
                    f.write(f"{location[0]} {location[1]} {location[2]} 0.3333333333 0.3333333333 0.000 {label}->K\n")
                    f.write(f"{location[0]} {location[1]} {location[2]} 0.500 0.000 0.000 {label}->M\n")
        else:
            print("Error: VBM location not found in the high symmetry path")


for mid in midlist:
    bandgappath = os.path.join(mid.strip(), 'band')
    vbmlocation, cbmlocation = get_vbm_cbm_locations(bandgappath, "")
    vbm_vpkitpath = os.path.join(mid.strip(), 'VPKIT.in')

    with open(vbm_vpkitpath, "w", newline="\n") as f:
        f.write("1\n6\n0.008\n2\n")

    writevpkit_kp(vbm_vpkitpath, vbmlocation, label="V")

    if np.linalg.norm([cbmlocation[i] - vbmlocation[i] for i in range(3)]) < 0.01:
        continue
    else:
        cbm_vpkitpath = os.path.join(mid.strip(), 'cbm_VPKIT.in')
        with open(cbm_vpkitpath, "w", newline="\n") as f:
            f.write("1\n6\n0.008\n2\n")
        writevpkit_kp(cbm_vpkitpath, cbmlocation, label="C")

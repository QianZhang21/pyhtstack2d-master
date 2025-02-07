"""For a more detailed analysis of PROCAR, it is necessary to refer to the code available at [
https://github.com/QijingZheng/pyband], which was developed by Associate Professor Qijing Zheng from the University
of Science and Technology of China. This code is solely used for analyzing band alignment and is thus a simplified
version. """
import os
import re
import numpy as np


def readPROCAR(propath, vbm_index=None, cbm_index=None, weightsave=False):
    if os.path.exists(os.path.join(propath, "weight.npy")) and os.path.exists(os.path.join(propath, "k_points.npy")):
        weights = np.load(os.path.join(propath, "weight.npy"))
        k_points = np.load(os.path.join(propath, "k_points.npy"))
    else:
        filepath = os.path.join(propath, "PROCAR")
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"The file {filepath} was not found!")

        with open(filepath, 'r') as f:
            file_lines = [line.strip() for line in f if line.strip()]
        nkpts, nbands, nions = map(int, re.findall(r'\d+', file_lines[1]))

        weights = []
        k_points = []
        for line in file_lines:
            if line.startswith("k-point"):
                match = re.search(r'k-point\s+\d+\s*:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+weight', line)
                if match:
                    coords = list(map(float, [match.group(1), match.group(2), match.group(3)]))
                    k_points.append(coords)
            if not re.search('[a-zA-Z]', line):
                weights.append(float(line.split()[-1]))

        k_points = np.array(k_points)
        weights = np.array(weights)

        nspin = weights.shape[0] // (nkpts * nbands * nions)
        assert nspin in [1, 2], f"Only support non-spin-polarized (1) or spin-polarized (2) calculations, " \
                                f"but got nspin={nspin}."

        weights = weights.reshape(nspin, nkpts, nbands, nions)
        if weightsave:
            np.save(os.path.join(os.path.dirname(filepath), 'weights.npy'), weights)
            np.save(os.path.join(os.path.dirname(filepath), 'k_points.npy'), k_points)

    if vbm_index is not None and cbm_index is not None:
        # vbmbandindex, cbmbandindex = vbm_index["band_index"][-1], cbm_index["band_index"][0]
        vbmbandindex = vbm_index["band_index"][-1] if type(vbm_index["band_index"]) is list else vbm_index["band_index"]
        cbmbandindex = cbm_index["band_index"][0] if type(cbm_index["band_index"]) is list else cbm_index["band_index"]
        spinindex = []
        closest_kpoint_index = []
        for index_dict in [vbm_index, cbm_index]:
            spinindex.append(0 if int(index_dict["spin"]) == 1 else 1)
            if "kpoint_coord" in index_dict:
                kpoint_coord = np.array(index_dict["kpoint_coord"])
                distances = np.linalg.norm(k_points - kpoint_coord, axis=1)
                closest_kpoint_index.append(np.argmin(distances))
            elif "kpoint_index" in index_dict:
                if type(index_dict["kpoint_index"]) is list:
                    closest_kpoint_index.append(index_dict["kpoint_index"][0])
                else:
                    closest_kpoint_index.append(index_dict["kpoint_index"])
            else:
                raise ValueError(f"Invalid kpointindex: {index_dict}")

        return weights[spinindex[0], closest_kpoint_index[0], vbmbandindex, :], weights[spinindex[1], closest_kpoint_index[1], cbmbandindex, :]
    else:
        return None, None

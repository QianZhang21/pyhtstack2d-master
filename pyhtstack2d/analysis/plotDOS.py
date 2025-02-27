import numpy as np
import os

import matplotlib as mpl

mpl.use('Agg')  # silent mode
from matplotlib import pyplot as plt

font = {'family': 'arial',
        'color': 'black',
        'weight': 'normal',
        'size': 22.0,
        }

color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
color_list2 = ['r', 'b', 'skyblue', 'steelblue', 'g', 'darkgreen', 'lightgreen', 'seagreen', 'turquoise', 'k']


class plotDOS(object):
    """
    The class is used to plot the projected Density-of-States of each element as calculated using VASP.
    path: str
        Path to the file containing the DOS calculation.
    e_range: tuple
        The energy range to plot the DOS data.
    skiptdos: bool
        Skip plotting the total DOS.
    """
    def __init__(self, path='./', e_range=(-3, 3)):
    # def __init__(self, path='./', e_range=(-3, 3), skiptdos=True):
        self.r_path = path
        self.E_min = e_range[0]
        self.E_max = e_range[1]
        # self.skiptdos = skiptdos
        self.plot_data_up = []
        self.plot_data_down = []
        self.plot_total_up = []
        self.plot_total_down = []
        self.labels = []
        self.read_dos_data()

    def read_dos_data(self):
        # if not self.skiptdos:
        #     if not any("TDOS.dat" in x for x in os.listdir(self.r_path)):
        #         os.system(f"(cd {self.r_path} && vaspkit -task 111 >/dev/null 2>&1)")
        #     if os.path.exists(os.path.join(self.r_path, "TDOS.dat")):
        #         tmp_data = np.array(np.loadtxt(os.path.join(self.r_path, "TDOS.dat")))
        #         self.plot_total_up.append(tmp_data[:, 1])
        #         self.plot_total_down.append(tmp_data[:, 2])
        #     else:
        #         self.skiptdos = True

        if not any("PDOS_" in x for x in os.listdir(self.r_path)):
            os.system(f"(cd {self.r_path} && vaspkit -task 113 >/dev/null 2>&1)")

        for filist in os.listdir(self.r_path):
            if "PDOS_" in filist and "IPDOS_" not in filist:
                tmp_data = np.array(np.loadtxt(os.path.join(self.r_path, filist)))
                self.plot_data_up.append(tmp_data[:, 0])
                self.plot_data_down.append(tmp_data[:, 0])
                break

        for filist in os.listdir(self.r_path):
            if "PDOS_" in filist:
                f_label = str(filist).split("_", 3)
                self.labels.append(f_label[1])

        self.labels = set(self.labels)
        for element in self.labels:
            for filist in os.listdir(self.r_path):
                if "IPDOS_" not in filist and "PDOS_" + element + "_UP" in filist:
                    tmp_data = np.array(np.loadtxt(os.path.join(self.r_path, filist)))
                    self.plot_data_up.append(tmp_data[:, -1])
                elif "IPDOS_" not in filist and "PDOS_" + element + "_DW" in filist:
                    tmp_data = np.array(np.loadtxt(os.path.join(self.r_path, filist)))
                    self.plot_data_down.append(tmp_data[:, -1])
                elif "IPDOS_" not in filist and "PDOS_" + element in filist:
                    tmp_data = np.array(np.loadtxt(os.path.join(self.r_path, filist)))
                    self.plot_data_up.append(tmp_data[:, -1])

    def plot_dos(self):
        plot_data_up = np.array(self.plot_data_up).T
        plot_data_dw = np.array(self.plot_data_down).T
        axe = plt.subplot(111)
        axe.axvline(x=0, linestyle='--', linewidth=1, color='0.75')
        emin_idy = next(x[0] for x in enumerate(plot_data_up[:, 0]) if x[1] >= self.E_min)
        emax_idy = len(plot_data_up[:, 0]) - next(
            y[0] for y in enumerate(reversed(plot_data_up[:, 0])) if y[1] <= self.E_max)
        if len(plot_data_dw[0]) > 1:
            dos_ymin = min(plot_data_dw[emin_idy:emax_idy + 1, 1] * 1.05)
        else:
            dos_ymin = 0
        dos_ymax = max([max(plot_data_up[emin_idy:emax_idy + 1, 1]) * 1.05, abs(dos_ymin)])
        for i, element in enumerate(self.labels):
            axe.plot(plot_data_up[:, 0], plot_data_up[:, i + 1], linewidth=1.0, label=str(element), color=color_list[i])
            if len(plot_data_dw[0]) > 1:
                axe.plot(plot_data_dw[:, 0], plot_data_dw[:, i + 1], linewidth=1.0, color=color_list[i])
                dos_ymin_tmp = min(plot_data_dw[emin_idy:emax_idy + 1, i + 1] * 1.05)
            else:
                dos_ymin_tmp = 0
            dos_ymax_tmp = max([max(plot_data_up[emin_idy:emax_idy + 1, i + 1]) * 1.05, abs(dos_ymin_tmp)])
            dos_ymin = min([dos_ymin, dos_ymin_tmp])
            dos_ymax = max([dos_ymax, dos_ymax_tmp])
        axe.set_xlabel(r'E-E$_\mathrm{f}$ (eV)', fontdict=font)
        axe.set_ylabel(r'DOS (a.u.)', fontdict=font)
        plt.legend()
        plt.xticks(fontsize=font['size'] - 2, fontname=font['family'])
        plt.yticks(fontsize=font['size'] - 2, fontname=font['family'])
        axe.set_xlim((self.E_min, self.E_max))
        axe.set_ylim((dos_ymin, dos_ymax))
        # axe.set_yticks([])
        fig = plt.gcf()
        fig.set_size_inches(8, 3)
        plt.savefig(os.path.join(self.r_path, 'dos.png'), dpi=300, bbox_inches='tight')

        plt.figure()
        axe_de = plt.subplot(111)
        axe_de.axhline(y=0, linestyle='--', linewidth=1, color='0.75')
        for i, element in enumerate(self.labels):
            axe_de.plot(plot_data_up[:, i + 1], plot_data_up[:, 0], linewidth=1.0, label=str(element),
                        color=color_list[i])
            if len(plot_data_dw[0]) > 1:
                axe_de.plot(plot_data_dw[:, i + 1], plot_data_dw[:, 0], linewidth=1.0, color=color_list[i])
        axe_de.set_xlabel(r'DOS (a.u.)', fontdict=font)
        axe_de.set_ylabel(r'E-E$_\mathrm{f}$ (eV)', fontdict=font)
        plt.legend()
        plt.xticks(fontsize=font['size'] - 2, fontname=font['family'])
        plt.yticks(fontsize=font['size'] - 2, fontname=font['family'])
        axe_de.set_ylim((self.E_min, self.E_max))
        axe_de.set_xlim((dos_ymin, dos_ymax))
        # axe_de.set_yticks([])
        fig_de = plt.gcf()
        fig_de.set_size_inches(3, 8)
        plt.savefig(os.path.join(self.r_path, 'dos_de.png'), dpi=300, bbox_inches='tight')

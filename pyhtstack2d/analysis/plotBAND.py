import os
import warnings

import matplotlib as mpl
import matplotlib.lines as mlines
import numpy as np

mpl.use('Agg')  # silent mode
from matplotlib.gridspec import GridSpec
from matplotlib import pyplot as plt
# plt.rc('font', family='Arial')  # set the default font family to Arial

color_list = ['g', 'c', 'm', 'y', 'k', 'Orange', 'Brown', 'Purple', 'Pink', 'Gray', 'Olive', 'Teal', 'Cyan',
              'Lime', 'Maroon', 'Navy', 'Indigo', 'DarkRed', 'DarkOrange', 'DarkGreen', 'DarkBlue', 'DarkCyan']


def run_vaspkit_task(bandpath, task, filename):
    """Helper function to run a vaspkit task and check for output file existence."""
    output_path = os.path.join(bandpath, filename)
    if not os.path.exists(output_path):
        os.system(f"(cd {bandpath} && vaspkit -task {task} >/dev/null 2>&1)")
    return os.path.exists(output_path)


class plotBS:
    def __init__(self, bandpath, savename=None, savepath=None, indices=None, natomlayer1=None, erange=(-3, 3),
                 elempro=True, hybrid=False, pmg=True, vaspkit=True,
                 BSDOSPlotter=None, BSPlotter=None, Vasprun=None, imag_format='pdf', dpi=300):
        """
        Plot band structure and density of states.

        Args:
            bandpath (str): Path to the band structure file.
            savename (str): Name of the saved plot.
            savepath (str): Path to save the plot.
            indices (list): Indices of the atoms in the layer.
            natomlayer1 (int): Number of atoms in the first layer.
            erange (tuple): Energy range for the plot.
            elempro (bool): Whether to plot the element projected band structure.
            hybrid (bool): Whether to plot the hybrid band structure.
            pmg (bool): Whether to use pymatgen to plot the band structure.
            vaspkit (bool): Whether to use vaspkit to generate the band.dat file.
            BSDOSPlotter (module): BSDOSPlotter module from pymatgen, when indices and natomlayer1 are None, this module
                will be used to plot the band structure and density of states.
            BSPlotter (module): BSPlotter module from pymatgen, when indices and natomlayer1 are provided, this module
                will be used to plot the layer projected band structure.
            Vasprun (module): Vasprun module from pymatgen.
            imag_format (str): Image format for the saved plot.
            dpi (int): Dots per inch for the saved plot.
        """
        self.ispin = False
        self.bandpath = bandpath
        if savename:
            savename = savename.rstrip("/")
            savename = savename.replace("/", "-")
        self.savename = f"{savename}-banddos.{imag_format}" if savename else f"banddos.{imag_format}"
        self.imag_format = imag_format
        self.dpi = dpi
        self.savepath = savepath if savepath else bandpath
        self.nlayer = 2 if natomlayer1 is not None and indices is not None else 1
        self.indices = indices
        self.natomlayer1 = natomlayer1
        self.erange = erange
        self.elempro = elempro
        self.hybrid = hybrid

        self.pmg = pmg
        self.vaspkit = vaspkit
        assert self.pmg or self.vaspkit, "Please choose either pymatgen or vaspkit to plot the band structure."
        assert os.path.exists(
            os.path.join(self.bandpath, 'vasprun.xml')), f"vasprun.xml not found in the directory {self.bandpath}."
        if self.pmg:
            assert (
                               BSDOSPlotter or BSPlotter) and Vasprun, "Please provide the BSDOSPlotter and Vasprun modules from pymatgen."
        self.BSDOSPlotter = BSDOSPlotter
        self.BSPlotter = BSPlotter
        self.Vasprun = Vasprun
        # if self.pmg:
        #     self.BSDOSPlotter = __import__('pymatgen.electronic_structure.plotter',
        #                                    fromlist=['BSDOSPlotter']).BSDOSPlotter
        #     self.Vasprun = __import__('pymatgen.io.vasp.outputs', fromlist=['Vasprun']).Vasprun
        # elif self.vaspkit:
        #     print("Please make sure that VASPKIT >= 1.3.1 is installed.")

        self.procar = False
        if os.path.exists(os.path.join(self.bandpath, 'PROCAR')):
            self.procar = True

    def plotbsdos(self):
        weight_l1, weight_l2, weight_dw_l1, weight_dw_l2 = None, None, None, None

        def plot_band_sactter(ax, band_data, weights, color, marker):
            for i in range(band_data.shape[1] - 1):
                ax.scatter(band_data[:, 0], band_data[:, 1 + i], s=weights[:, i] * 60,
                           marker=marker, color=color, alpha=0.5, zorder=1, lw=0.0)

        if self.pmg:
            force_hybrid_mode = True if self.hybrid else False
            bs_vasprun = self.Vasprun(os.path.join(self.bandpath, 'vasprun.xml'), parse_projected_eigen=True)
            bs_data = bs_vasprun.get_band_structure(line_mode=True, force_hybrid_mode=force_hybrid_mode)
            if self.nlayer == 2:
                projections = bs_data.projections
                weight_l1, weight_l2, weight_dw_l1, weight_dw_l2 = self.projeclayer(projections)
                legendist = []
                if self.BSPlotter:
                    bsplot = self.BSPlotter(bs_data)
                    bsplot.get_plot()
                    bsplotdata = bsplot.bs_plot_data()
                    kpoint = bsplotdata['distances']
                    band_energies = bsplotdata['energy']['1']
                    band_energies_dw = bsplotdata['energy']['-1'] if '-1' in bsplotdata['energy'].keys() else None
                    kpoint = [np.array(segment) for segment in kpoint]
                    banddata = [np.array(segment) for segment in band_energies]
                    banddata_dw = [np.array(segment) for segment in band_energies_dw] if band_energies_dw is not None else None
                    kpoint = np.concatenate(kpoint, axis=0)
                    banddata = np.concatenate(banddata, axis=1).T

                    def legendlabel(color, marker, markersize, linestyle, label):
                        return mlines.Line2D([], [], color=color, marker=marker, markersize=markersize, linestyle=linestyle, label=label)

                    if weight_l1 is not None:
                        for i in range(banddata.shape[1]):
                            plt.scatter(kpoint, banddata[:, i], s=weight_l1[:, i] * 150,
                                        marker='o', color='red', alpha=0.5, zorder=1, lw=0.0)
                            plt.scatter(kpoint, banddata[:, i], s=weight_l2[:, i] * 150,
                                        marker='o', color='blue', alpha=0.5, zorder=1, lw=0.0)
                        if banddata_dw is not None:
                            banddata_dw = np.concatenate(banddata_dw, axis=0)
                            for i in range(banddata_dw.shape[1]):
                                plt.scatter(kpoint, banddata_dw[:, i], s=weight_dw_l1[:, i] * 150,
                                            marker='*', color='red', alpha=0.5, zorder=1, lw=0.0)
                                plt.scatter(kpoint, banddata_dw[:, i], s=weight_dw_l2[:, i] * 150,
                                            marker='*', color='blue', alpha=0.5, zorder=1, lw=0.0)
                            legendist.append(legendlabel('red', 'o', 12, 'None', 'L1-UP'))
                            legendist.append(legendlabel('blue', 'o', 12, 'None', 'L2-UP'))
                            legendist.append(legendlabel('red', '*', 12, 'None', 'L1-DW'))
                            legendist.append(legendlabel('blue', '*', 12, 'None', 'L2-DW'))
                        else:
                            legendist.append(legendlabel('red', 'o', 12, 'None', 'L1'))
                            legendist.append(legendlabel('blue', 'o', 12, 'None', 'L2'))
                    if weight_l1 is not None or banddata_dw is not None:
                        plt.legend(fontsize=24, handles=legendist, loc='best')
                    plt.savefig(os.path.join(self.savepath, self.savename), format=self.imag_format, dpi=self.dpi)
                    plt.close()
                else:
                    warnings.warn(f"Band structure plotter is not available. "
                                  f"Please install the package 'pymatgen' and import 'BSPlotter' or 'BSDOSPlotter'.")
            else:
                dos_vasprun = self.Vasprun(os.path.join(self.bandpath, 'vasprun.xml'))
                dos_data = dos_vasprun.complete_dos

                bspro = "elements" if self.procar and self.elempro else None
                dospro = "elements" if self.procar and self.elempro else None
                banddos_fig = self.BSDOSPlotter(bs_projection=bspro, dos_projection=dospro,
                                                vb_energy_range=abs(self.erange[0]),
                                                cb_energy_range=abs(self.erange[1]))
                banddos_fig.get_plot(bs=bs_data, dos=dos_data)
                # banddos_fig.save_plot(os.path.join(self.savepath, self.savename))
                plt.savefig(os.path.join(self.savepath, self.savename), format=self.imag_format, dpi=self.dpi)
                plt.close()

        elif self.vaspkit:
            if self.nlayer == 2:
                weight_l1, weight_l2, weight_dw_l1, weight_dw_l2 = self.projeclayer()

            suffix = ["_UP", "_DW"]
            banddata_dw = None
            taskcode = 252 if self.hybrid else 211
            if not run_vaspkit_task(self.bandpath, taskcode, "REFORMATTED_BAND.dat"):
                self.ispin = True
                for suf in suffix:
                    if not os.path.exists(os.path.join(self.bandpath, f"REFORMATTED_BAND{suf}.dat")):
                        raise FileNotFoundError(f"Band data file not found in the directory {self.bandpath}.")
                banddata = np.loadtxt(os.path.join(self.bandpath, f"REFORMATTED_BAND_UP.dat"), dtype=np.float64)
                banddata_dw = np.loadtxt(os.path.join(self.bandpath, f"REFORMATTED_BAND_DW.dat"), dtype=np.float64)
            else:
                banddata = np.loadtxt(os.path.join(self.bandpath, "REFORMATTED_BAND.dat"), dtype=np.float64)

            klabels = os.path.join(self.bandpath, "KLABELS")
            group_labels, xtick = [], []
            if os.path.exists(klabels):
                with open(klabels, 'r', encoding='utf-8') as f:
                    lines = f.readlines()[1:]
                for line in lines:
                    tokens = line.split()
                    if len(tokens) == 2 and not line.startswith('*'):
                        group_labels.append(tokens[0])
                        xtick.append(float(tokens[1]))
                group_labels = [u'Γ' if label == 'GAMMA' else label for label in group_labels]

            tdosdata = None
            if self.procar:
                tdosdata_exists = run_vaspkit_task(self.bandpath, 111, "TDOS.dat")
                if tdosdata_exists:
                    tdosdata = np.loadtxt(os.path.join(self.bandpath, "TDOS.dat"), dtype=np.float64)

            if tdosdata is None:
                fig, bs_ax = plt.subplots()
            else:
                gs = GridSpec(1, 2, width_ratios=[2, 1])
                bs_ax, dos_ax = plt.subplot(gs[0]), plt.subplot(gs[1])

            bs_ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1, color='0.75')
            if len(xtick) != 0:
                bs_ax.set_xticks(xtick)
                bs_ax.set_xlim((xtick[0], xtick[-1]))
                bs_ax.set_xticklabels(group_labels, rotation=0)
                for i in xtick[1:-1]:
                    bs_ax.axvline(x=i, ymin=0, ymax=1, linestyle='--', linewidth=1, color='0.75')

            if weight_l1 is not None:
                plot_band_sactter(bs_ax, banddata, weight_l1, 'red', 'o')
                plot_band_sactter(bs_ax, banddata, weight_l2, 'blue', 'o')
            else:
                bs_ax.plot(banddata[:, 0], banddata[:, 1], linewidth=2, color='black', label='Spin Up')
                bs_ax.plot(banddata[:, 0], banddata[:, 2:], linewidth=2, color='black')

            if banddata_dw is not None:
                if weight_dw_l1 is not None:
                    plot_band_sactter(bs_ax, banddata_dw, weight_dw_l1, 'red', '*')
                    plot_band_sactter(bs_ax, banddata_dw, weight_dw_l2, 'blue', '*')
                    bs_ax.plot([], [], 'o', color='red', label='L1-UP')
                    bs_ax.plot([], [], 'o', color='blue', label='L2-UP')
                    bs_ax.plot([], [], '*', color='red', label='L1-DW')
                    bs_ax.plot([], [], '*', color='blue', label='L2-DW')
                else:
                    bs_ax.plot(banddata_dw[:, 0], banddata_dw[:, 1], linewidth=2, color='red', label='Spin Down')
                    bs_ax.plot(banddata_dw[:, 0], banddata_dw[:, 2:], linewidth=2, color='red')
            else:
                bs_ax.plot([], [], 'o', color='red', label='L1')
                bs_ax.plot([], [], 'o', color='blue', label='L2')

            bs_ax.set_ylim((self.erange[0], self.erange[1]))
            # bs_ax.tick_params(direction='in')
            bs_ax.set_ylabel(r'E-E$_\mathrm{f}$ (eV)', fontsize=14)
            yticks = np.arange(self.erange[0], self.erange[1] + 0.5, 0.5)
            bs_ax.set_yticks(yticks)
            bs_ax.set_yticklabels([f"{tick:.1f}" for tick in yticks])
            # bs_ax.set_yticklabels(np.arange(self.erange[0], self.erange[1] + 1e-8, 0.5))
            if self.ispin or weight_l1 is not None:
                bs_ax.legend(loc='lower right')

            if tdosdata is not None:
                dos_energy = tdosdata[:, 0]
                dos = tdosdata[:, 1]
                dos_ax.plot(dos, dos_energy, color=(0.6, 0.6, 0.6), linewidth=2, label='Total')
                dos_ax.fill_betweenx(dos_energy, 0, dos, where=(dos > 0), color=(0.7, 0.7, 0.7),
                                     facecolor=(0.7, 0.7, 0.7))
                if tdosdata.shape[1] == 3:
                    dos_down = tdosdata[:, 2]
                    dos_ax.plot(dos_down, dos_energy, color=(0.6, 0.6, 0.6), linewidth=2)
                    dos_ax.fill_betweenx(dos_energy, dos_down, 0, color=(0.7, 0.7, 0.7), facecolor=(0.7, 0.7, 0.7))

                elements, pdosdata, pdosdata_dw = self.readPDOS()
                if pdosdata is not None:
                    for i, elem_ in enumerate(elements):
                        dos_ax.plot(pdosdata[i], dos_energy, linewidth=1.5, label=elem_, color=color_list[i])
                        if pdosdata_dw is not None:
                            dos_ax.plot(pdosdata_dw[i], dos_energy, linewidth=1.5, color=color_list[i])

                emin_idx = next(x[0] for x in enumerate(dos_energy) if x[1] >= self.erange[0])
                emax_idx = len(dos_energy) - next(
                    x[0] for x in enumerate(reversed(dos_energy)) if x[1] <= self.erange[1])
                dos_xmin = (0 if tdosdata.shape[1] == 2 else min(dos_down[emin_idx: emax_idx + 1] * 1.05))
                dos_xmax = max([max(dos[emin_idx:emax_idx]) * 1.05, abs(dos_xmin)])
                dos_ax.set_xlim(dos_xmin, dos_xmax)
                dos_ax.set_xticklabels([])
                dos_ax.axhline(y=0, xmin=0, xmax=1, linestyle='--', linewidth=1, color='0.75')
                # dos_ax.hlines(y=0, xmin=dos_xmin, xmax=dos_xmax, linestyle='--', linewidth=1, color='0.75')
                dos_ax.set_xlabel("DOS", fontsize=14)
                dos_ax.legend(fancybox=True, loc='upper right')  # dos_ax.legend(loc=(0.4, 0.2))
                dos_ax.set_ylim((self.erange[0], self.erange[1]))
                dos_ax.set_yticks([])
                dos_ax.set_yticklabels([])
                dos_ax.set_yticks(np.arange(self.erange[0], self.erange[1] + 1e-5, 1))
                dos_ax.grid(color=[0.75, 0.75, 0.75], linestyle="--", linewidth=1)
                # dos_ax.tick_params(direction='in')
                plt.subplots_adjust(wspace=0.1)
                # fig = plt.gcf()
                # fig.set_size_inches((9, 8.5))
                plt.savefig(os.path.join(self.savepath, self.savename), bbox_inches='tight', format=self.imag_format, dpi=self.dpi)
                plt.close()
            else:
                plt.subplots_adjust(wspace=0.1)
                plt.savefig(os.path.join(self.savepath, self.savename), bbox_inches='tight', format=self.imag_format,
                            dpi=self.dpi)
                plt.close()

    def readPDOS(self):
        if not self.elempro or not self.procar:
            return None, None
        os.system(f"(cd {self.bandpath} && vaspkit -task 113 >/dev/null 2>&1)")
        elements = []
        for file_ in os.listdir(self.bandpath):
            if file_.startswith("PDOS_"):
                elements.append(file_.split("_")[1].split(".")[0])
        pdosdata, pdosdata_dw = [], []
        elements = list(set(elements))
        for elem_ in elements:
            if self.ispin:
                pdosdata.extend([np.loadtxt(os.path.join(self.bandpath, f"PDOS_{elem_}_UP.dat"))[:, -1]])
                pdosdata_dw.extend([np.loadtxt(os.path.join(self.bandpath, f"PDOS_{elem_}_DW.dat"))[:, -1]])
            else:
                pdosdata.extend([np.loadtxt(os.path.join(self.bandpath, f"PDOS_{elem_}.dat"))[:, -1]])
        pdosdata = np.array(pdosdata)
        pdosdata_dw = np.array(pdosdata_dw) if len(pdosdata_dw) > 0 else None
        return elements, pdosdata, pdosdata_dw

    def projeclayer(self, bs_projections=None):
        weightfile = os.path.join(self.bandpath, 'weights.npy')
        atomlayer1 = [index for index, i_ in enumerate(self.indices) if i_ < self.natomlayer1]
        atomlayer2 = [index for index, i_ in enumerate(self.indices) if i_ >= self.natomlayer1]

        if os.path.exists(weightfile):
            weightread = np.load(weightfile)
            weightread_l1 = np.sum(weightread[:, :, :, atomlayer1], axis=3)
            weightread_l2 = np.sum(weightread[:, :, :, atomlayer2], axis=3)
            weight_l1 = weightread_l1[0]
            weight_l2 = weightread_l2[0]
            if weightread.shape[0] == 2:
                weight_dw_l1 = weightread_l1[1]
                weight_dw_l2 = weightread_l2[1]
            else:
                weight_dw_l1, weight_dw_l2 = None, None
            return weight_l1, weight_l2, weight_dw_l1, weight_dw_l2
        else:
            if self.pmg and bs_projections is not None:
                weight_tmp = {}
                key_tmp = []
                for spin_key, spin_projections in bs_projections.items():
                    key_tmp.append(str(spin_key))
                    weight_tmp[str(spin_key)] = np.sum(spin_projections, axis=2)
                weight_l1 = np.sum(weight_tmp[key_tmp[0]][:, :, atomlayer1], axis=2).T
                weight_l2 = np.sum(weight_tmp[key_tmp[0]][:, :, atomlayer2], axis=2).T
                if len(key_tmp) == 2:
                    weight_dw_l1 = np.sum(weight_tmp[key_tmp[1]][:, :, atomlayer1], axis=2).T
                    weight_dw_l2 = np.sum(weight_tmp[key_tmp[1]][:, :, atomlayer2], axis=2).T
                else:
                    weight_dw_l1, weight_dw_l2 = None, None
                return weight_l1, weight_l2, weight_dw_l1, weight_dw_l2
            elif self.vaspkit:
                atomlayer1 = " ".join([str(i_ + 1) for i_ in atomlayer1])
                atomlayer2 = " ".join([str(i_ + 1) for i_ in atomlayer2])
                os.system(f"echo -e '214\\n1\\n{atomlayer1}\\n' | (cd {self.bandpath} && vaspkit >/dev/null 2>&1)")

                def parse_dat_file(filepath):
                    nkpoints, nbands = None, None
                    with open(filepath, 'r') as file:
                        for line in file:
                            if line.startswith("# NKPTS & NBANDS:"):
                                parts = line.split(':')[-1].split()
                                nkpoints = int(parts[0])
                                nbands = int(parts[1])
                                break
                    return nkpoints, nbands

                def organize_weight(weightfile, nkpoints, nbands):
                    weight = np.loadtxt(weightfile)[:, -1].reshape(nbands, nkpoints)
                    for i in range(nbands):
                        if i % 2 != 0:
                            weight[i] = weight[i][::-1]
                    return weight

                pband_sum_path = os.path.join(self.bandpath, 'PBAND_SUM.dat')
                if os.path.exists(pband_sum_path):
                    nkpoints, nbands = parse_dat_file(pband_sum_path)
                    weight_l1 = organize_weight(pband_sum_path, nkpoints, nbands).T
                    # weight_l1 = np.loadtxt(pband_sum_path)[:, -1].reshape(nbands, nkpoints)
                    # for i in range(nbands):
                    #     if i % 2 != 0:
                    #         weight_l1[i] = weight_l1[i][::-1]
                    os.remove(pband_sum_path)
                    os.system(f"echo -e '214\\n1\\n{atomlayer2}\\n' | (cd {self.bandpath} && vaspkit >/dev/null 2>&1)")
                    if os.path.exists(pband_sum_path):
                        # weight_l2 = np.loadtxt(pband_sum_path)[:, -1].reshape(nbands, nkpoints).T
                        weight_l2 = organize_weight(pband_sum_path, nkpoints, nbands).T
                        os.remove(pband_sum_path)
                        return weight_l1, weight_l2, None, None

                pband_sum_up_path = os.path.join(self.bandpath, 'PBAND_SUM_UP.dat')
                pband_sum_dw_path = os.path.join(self.bandpath, 'PBAND_SUM_DW.dat')
                if os.path.exists(pband_sum_up_path):
                    nkpoints, nbands = parse_dat_file(pband_sum_up_path)
                    # weight_l1 = np.loadtxt(pband_sum_up_path)[:, -1].reshape(nbands, nkpoints).T
                    # weight_dw_l1 = np.loadtxt(pband_sum_dw_path)[:, -1].reshape(nbands, nkpoints).T
                    weight_l1 = organize_weight(pband_sum_up_path, nkpoints, nbands).T
                    weight_dw_l1 = organize_weight(pband_sum_dw_path, nkpoints, nbands).T
                    os.remove(pband_sum_up_path)
                    os.remove(pband_sum_dw_path)
                    os.system(f"echo -e '214\\n1\\n{atomlayer2}\\n' | (cd {self.bandpath} && vaspkit >/dev/null 2>&1)")
                    if os.path.exists(pband_sum_up_path):
                        # weight_l2 = np.loadtxt(pband_sum_up_path)[:, -1].reshape(nbands, nkpoints).T
                        # weight_dw_l2 = np.loadtxt(pband_sum_dw_path)[:, -1].reshape(nbands, nkpoints).T
                        weight_l2 = organize_weight(pband_sum_up_path, nkpoints, nbands).T
                        weight_dw_l2 = organize_weight(pband_sum_dw_path, nkpoints, nbands).T
                        os.remove(pband_sum_up_path)
                        os.remove(pband_sum_dw_path)
                        return weight_l1, weight_l2, weight_dw_l1, weight_dw_l2
                return None, None, None, None

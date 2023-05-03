# -*- coding: utf-8 -*-
"""
Tools to analyze the admittance matrix data
"""

import plotting as plt
import numpy as np
import scipy.interpolate
import re as regex

class Admittance(object):
    """
    Class for handling the admittance data
    """
    def __init__(self, *args, **kwargs):

        # Default class conditions
        self.all_data_loaded = False
        self.units = r'$C$ [fF]'
        self.scale = 1e-15
        self.fscale = 1e9
        self.symmetrize = False

        # Update the arguments and keyword arguments
        self.__dict__.update(locals())
        for k, v in kwargs.items():
            setattr(self, k, v)


    def __del__(self):
        pass

    def load_all_data_from_file(self, fname, skip_header=5):
        """
        Loads the admittance data from a single file using the header
        data to extract the gate voltages, frequencies
        """
        # Read the data and order into Vg-labeled dictionary
        data       = np.genfromtxt(fname, dtype=complex,
                                   skip_header=skip_header).T
        # Zero out the nans
        data[np.isnan(data)] = 0
        self.freq  = np.unique(data[4].real)
        Nf         = self.freq.size
        Vg_uniq    = np.unique(data[0].real)
        self.Vg    = Vg_uniq
        Nvg        = Vg_uniq.size
        Y          = data[6:]
        self.Nt = int(np.sqrt(Y.shape[0]))
        self.data = {}
        Y = np.reshape(Y, [self.Nt * self.Nt, Nvg, self.Nt, Nf])
        for vidx, vg in enumerate(Vg_uniq):
            for fidx, f in enumerate(self.freq):
                Y_key = f'y_vg_{vg}_f_{f/self.fscale}GHz'.replace('.', 'p')
                Yin = np.zeros(self.Nt * self.Nt, dtype=Y.dtype)
                for i in range(self.Nt):
                    # Yin += Y[vidx, i, fidx, :]
                    # Yin += Y[:, fidx, i, vidx]
                    Yin += Y[:, vidx, i, fidx]
                self.data[Y_key] = [vg, Yin.reshape([self.Nt, self.Nt])]

        self.all_data_loaded = True

    def group_by_frequency_at_vg(self, Vg):
        """
        Group the admittance data at a particular Vg by frequency
        """
        try:
            Y = np.array([self.data[f'y_vg_{Vg}_f_{f/self.fscale}GHz'\
                .replace('.', 'p')][1] for f in self.freq])
        except KeyError as err:
            print(f'{err}\nValid keys: {self.data.keys()}')

        return Y

    def plot_yij_vs_f_at_vg(self, i, j, Vg, fname):
        """
        Plot the admittance matrix elements as a function of gate voltage
        """
        # Read the data into an array
        Y   = self.group_by_frequency_at_vg(Vg)
        Yij = np.array([YY[i, j] for YY in Y])

        # Initialize the figure and axes
        myplt = plt.MPLPlotWrapper()
        myplt.plot(self.freq / self.fscale, np.real(Yij), 'o-',
                   color='tab:blue')
        ax2 = myplt.ax.twinx()
        ax2.plot(self.freq / self.fscale, np.imag(Yij), 's-',
                 color='tab:red')
        myplt.xlabel = r'$f$ [GHz]'
        myplt.ylabel = r'$\mathrm{Re}Y_{%d%d}$ [S]' % (i + 1, j + 1)
        ax2.set_ylabel(r'$\mathrm{Im}Y_{%d%d}$ [S]' % (i + 1, j + 1),
                fontsize=myplt.fsize)
        myplt.set_axes_fonts(ax=ax2)

        # Write the figure to file
        myplt.write_fig_to_file(fname)

    def plot_yij_vs_f_vs_vg(self, i, j, fname, real_imag='real'):
        """
        Plot the admittance matrix elements as a function of gate voltage
        """
        # Initialize the figure
        myplt = plt.MPLPlotWrapper()
        mkrs = myplt.marker_cycle
        mlen = len(mkrs)

        # Read the data into an array
        for vidx, Vg in enumerate(self.Vg):

            # Get the data at Vg
            Y   = self.group_by_frequency_at_vg(Vg)
            Yij = np.array([YY[i, j] for YY in Y])

            # Initialize the figure and axes
            if real_imag == 'real':
                myplt.plot(self.freq / self.fscale, np.real(Yij),
                            marker=mkrs[vidx % mlen],
                            label=r'$V_g=$%.1f V' % Vg)
                myplt.ylabel = r'$\mathrm{Re}Y_{%d%d}$ [S]' % (i + 1, j + 1)
            elif real_imag == 'imag':
                myplt.plot(self.freq / self.fscale, np.imag(Yij), 
                            marker=mkrs[vidx % mlen],
                            label=r'$V_g=$%.1f V' % Vg)
                myplt.ylabel = r'$\mathrm{Im}Y_{%d%d}$ [S]' % (i + 1, j + 1)

        myplt.xlabel = r'$f$ [GHz]'

        # Write the figure to file
        myplt.set_leg_outside()
        myplt.write_fig_to_file(fname)

    def plot_y_all_vs_vg(self, fname):
        """
        Plot the full capacitance matrix as a function of gate voltage
        """
        # Read the data into an array
        Yij = np.array([Yval[1] for Yval in self.data.values()])
        Vg  = np.array([Yval[0] for Yval in self.data.values()])

        # Initialize the figure and axes
        myplt = plt.MPLPlotWrapper()
        idx = 0
        mrks = myplt.marker_cycle
        mlen = len(mrks)
        for i in range(self.Nt):
            for j in range(i, self.Nt):
                print(f'Y[{i}, {j}]')
                myplt.plot(Vg, np.abs(Yij[:, i, j]),
                           marker=mrks[idx % mlen], ls='-',
                    label=r'$Y_{%d%d}$' % (i + 1, j + 1))
                idx += 1

        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$|Y_{ij}|$ [fF]'
        myplt.set_leg_outside()

        # Write the figure to file
        myplt.write_fig_to_file(fname)

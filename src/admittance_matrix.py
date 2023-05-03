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
        Nvg        = Vg_uniq.size
        Y          = data[6:]
        self.Nterm = int(np.sqrt(Y.shape[0]))
        self.data = {}
        for vidx, vg in enumerate(Vg_uniq):
            for fidx, f in enumerate(self.freq):
                Y_key = f'ymatrix_vg_{vg}_f_{f/self.fscale}GHz'.replace('.', 'p')
                Yavg = 0
                for i in range(self.Nterm):
                    idx = vidx + Nvg * fidx + Nvg * Nf * i
                    iidd = vidx*Nf*self.Nterm + i
                    Yavg +=  Y[:, idx]
                Yin = Yavg.reshape([self.Nterm, self.Nterm]) # / self.Nterm
                if self.symmetrize:
                    Yin = 0.5 * (Yin + Yin.T)
                self.data[Y_key] = [vg, Yin]

        self.all_data_loaded = True

    def plot_yij_vs_vg(self, i, j, fname):
        """
        Plot the admittance matrix elements as a function of gate voltage
        """
        # Read the data into an array
        Yij = np.array([Yval[1][i, j] for Yval in self.data.values()])
        Vg  = np.array([Yval[0] for Yval in self.data.values()])

        # Initialize the figure and axes
        myplt = plt.MPLPlotWrapper()
        print(f'Y[{i}, {j}]')
        myplt.plot(Vg, np.abs(Yij), 'o-')
        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$|Y_{%d%d}|$ [S]' % (i + 1, j + 1)

        # Write the figure to file
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
        for i in range(self.Nterm):
            for j in range(i, self.Nterm):
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

# -*- coding: utf-8 -*-
"""
Tools to analyze the Maxwell capacitance matrix data
"""

import plotting as plt
import numpy as np
import scipy.interpolate
import re as regex

class MaxwellCapacitance(object):
    """
    Class for handling the Maxwell Capacitance data
    """
    def __init__(self, *args, **kwargs):

        # Default class conditions
        self.all_data_loaded = False
        self.units = r'$C$ [fF]'
        self.scale = 1e-15
        self.symmetrize = False

        # Update the arguments and keyword arguments
        self.__dict__.update(locals())
        for k, v in kwargs.items():
            setattr(self, k, v)


    def __del__(self):
        pass

    def load_all_data_from_file(self, fname, skip_header=5):
        """
        Loads the charge concentration data from a single file using the header
        data to extract the gate voltages
        """
        # Read the data and order into Vg-labeled dictionary
        data = np.genfromtxt(fname, skip_header=skip_header).T
        Vg        = data[0]
        Vg_uniq   = np.unique(Vg)[::-1]
        self.Nt   = Vg.size // Vg_uniq.size
        Nvg       = Vg_uniq.size
        C         = data[5:] / self.scale
        self.data = {}
        C         = np.reshape(C, [self.Nt * self.Nt, Nvg, self.Nt])
        for vidx, vg in enumerate(Vg_uniq):
            C_key = f'cmatrix_vg_{vg}'.replace('.', 'p')
            Cavg = 0
            for i in range(self.Nt):
                # Cavg +=  C[:, vidx * self.Nt + i]
                Cavg +=  C[:, vidx, i]
            Cin = Cavg.reshape([self.Nt, self.Nt]) # / self.Nt
            if self.symmetrize:
                Cin = 0.5 * (Cin + Cin.T)
            self.data[C_key] = [vg, Cin]

        self.all_data_loaded = True

    def plot_cij_vs_vg(self, i, j, fname):
        """
        Plot the capacitance matrix elements as a function of gate voltage
        """
        # Read the data into an array
        Cij = np.array([Cval[1][i, j] for Cval in self.data.values()])
        Vg  = np.array([Cval[0] for Cval in self.data.values()])

        # Initialize the figure and axes
        myplt = plt.MPLPlotWrapper()
        print(f'C[{i}, {j}]')
        myplt.plot(Vg, np.abs(Cij), 'o-')
        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$|C_{%d%d}|$ [fF]' % (i + 1, j + 1)

        # Write the figure to file
        myplt.write_fig_to_file(fname)

    def plot_c_all_vs_vg(self, fname, ylim=None, xrot=None):
        """
        Plot the full capacitance matrix as a function of gate voltage
        """
        # Read the data into an array
        Cij = np.array([Cval[1] for Cval in self.data.values()])
        Vg  = np.array([Cval[0] for Cval in self.data.values()])

        # Initialize the figure and axes
        myplt = plt.MPLPlotWrapper()
        idx = 0
        mrks = myplt.marker_cycle
        mlen = len(mrks)
        for i in range(self.Nt):
            # for j in range(i, self.Nt):
            for j in range(self.Nt):
                print(f'C[{i}, {j}]')
                myplt.plot(Vg, Cij[:, i, j],
                           marker=mrks[idx % mlen], ls='-',
                    label=r'$C_{%d%d}$' % (i + 1, j + 1))
                idx += 1

        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$C_{ij}$ [fF]'
        if xrot is not None:
            myplt.set_xaxis_rot(xrot)
        if ylim is not None:
            myplt.ylim = ylim
        myplt.set_leg_outside()

        # Write the figure to file
        myplt.write_fig_to_file(fname)

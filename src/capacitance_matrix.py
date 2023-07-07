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
        self.escale = 1e6
        self.C1 = 96.6e-15
        self.C2 = 96.6e-15

        # Update the arguments and keyword arguments
        self.__dict__.update(locals())
        for k, v in kwargs.items():
            setattr(self, k, v)


    def __del__(self):
        pass

    def load_all_data_from_file(self, fname, skip_header=5, vg_start_idx=None):
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
        self.data = {}
        C         = data[5:]
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

    def get_coupling_matrix(self, parasitic=False, inv=True):
        """
        Computes the coupling capacitance matrix
        """
        # Get the full C-matrix vs. Vg
        Cij = np.array([Cval[1] for Cval in self.data.values()])
        Vg  = np.array([Cval[0] for Cval in self.data.values()])

        if parasitic:
            # Compute the parasitic capacitance matrix
            C11 = self.C1 + Cij[:, 0, 0] + abs(Cij[:, 0, 1]) + abs(Cij[:, 0, 2])
            C12 = -abs(Cij[:, 0, 1]) 
            C13 = -abs(Cij[:, 0, 2])
            C22 = self.C2 + Cij[:, 1, 1] + abs(Cij[:, 0, 1]) + abs(Cij[:, 1, 2])
            C23 = -abs(Cij[:, 1, 2])
            C33 = Cij[:, 2, 2] + abs(Cij[:, 0, 2]) + abs(Cij[:, 1, 2])
            Cfull = np.array([[[C11[i], C12[i], C13[i]],
                               [C12[i], C22[i], C23[i]],
                               [C13[i], C23[i], C33[i]]] \
                                       for i in range(len(Vg))])
        else:
            Cfull = np.array([
                [[self.C1 + abs(Cij[i, 0, 1]), -abs(Cij[i, 0, 1])],
                 [-abs(Cij[i, 0, 1]), self.C2 + abs(Cij[i, 0, 1])]]
                for i in range(len(Vg))])

        # Return the original matrix or the inverse
        if inv:
            return np.array([np.linalg.inv(Cfull[i]) for i in range(len(Vg))])
        else:
            return Cfull

    def plot_cij_parasitic_simplified_compare(self, i, j, fname, 
                                              xrot=None, ylim=None,
                                              vg_start_idx=None):
        """
        Plot the Cij^-1 matrix elements, parasitic vs. simplified
        """
        # Physical constants to convert from mks to MHz
        e = 1.602176634e-19
        h = 6.62607015e-34

        # Get the full C-matrix vs. Vg
        Vg  = np.array([Cval[0] for Cval in self.data.values()])

        # Common plotting functions to both
        myplt = plt.MPLPlotWrapper()
        idx = 0
        mrks = myplt.marker_cycle
        mlen = len(mrks)

        # Compute the coupling matrices
        Cinvfull = self.get_coupling_matrix(parasitic=True)
        Cinvsimp = self.get_coupling_matrix(parasitic=False)
        Cfull = self.get_coupling_matrix(parasitic=True, inv=False)
        Csimp = self.get_coupling_matrix(parasitic=False, inv=False)

        # Plot the coupling matrices
        if vg_start_idx is not None:
            myplt.plot(Vg[vg_start_idx:],
              0.5 * e**2 * Cinvsimp[vg_start_idx:, i, j] / h / self.escale,
                       marker=mrks[0], ls='-',
                label=r'$C^{-1}_{\mathrm{simp}, %d%d}$' % (i + 1, j + 1))
            myplt.plot(Vg[vg_start_idx:],
              0.5 * e**2 * Cinvfull[vg_start_idx:, i, j] / h / self.escale,
                       marker=mrks[1], ls='-',
                label=r'$C^{-1}_{\mathrm{para}, %d%d}$' % (i + 1, j + 1))
        else:
            myplt.plot(Vg, 0.5 * e**2 * Cinvsimp[:, i, j] / h / self.escale,
                       marker=mrks[0], ls='-',
                label=r'$C^{-1}_{\mathrm{simp}, %d%d}$' % (i + 1, j + 1))
            myplt.plot(Vg, 0.5 * e**2 * Cinvfull[:, i, j] / h / self.escale,
                       marker=mrks[1], ls='-',
                label=r'$C^{-1}_{\mathrm{para}, %d%d}$' % (i + 1, j + 1))

        # Final plotting command common to both
        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$\frac{1}{2}e^2C^{-1}_{ij}$ [MHz]'
        if xrot is not None:
            myplt.set_xaxis_rot(xrot)
        if ylim is not None:
            myplt.ylim = ylim
        myplt.set_leg_outside()

        # Compute the on/off ratios
        if vg_start_idx is not None:
            Cinvfullmax   = Cinvfull[vg_start_idx, :, :]
            Cinvfullmin   = Cinvfull[-1, :, :]
            Cinvfullonoff = Cinvfullmax[i, j] / Cinvfullmin[i, j]
            Cinvsimpmax   = Cinvsimp[vg_start_idx, :, :]
            Cinvsimpmin   = Cinvsimp[-1, :, :]
            Cinvsimponoff = Cinvsimpmax[i, j] / Cinvsimpmin[i, j]

            Cij = np.array([Cval[1] for Cval in self.data.values()])
            Cijmax = Cij[vg_start_idx, :, :]
            Cijmin = Cij[-1, :, :]
            Cijonoff = Cijmax[i, j] / Cijmin[i, j]

        else:
            Cinvfullmax   = Cinvfull[0, :, :]
            Cinvfullmin   = Cinvfull[-1, :, :]
            Cinvfullonoff = Cinvfullmax[i, j] / Cinvfullmin[i, j]
            Cinvsimpmax   = Cinvsimp[0, :, :]
            Cinvsimpmin   = Cinvsimp[-1, :, :]
            Cinvsimponoff = Cinvsimpmax[i, j] / Cinvsimpmin[i, j]

            Cij = np.array([Cval[1] for Cval in self.data.values()])
            Cijmax = Cij[0, :, :]
            Cijmin = Cij[-1, :, :]
            Cijonoff = Cijmax[i, j] / Cijmin[i, j]

        np.set_printoptions(precision=2)
        print(f'Cfullmax^-1:\n{0.5 * e**2 * Cinvfullmax / h / self.escale}')
        print(f'Cfullmin^-1:\n{0.5 * e**2 * Cinvfullmin / h / self.escale}')
        print(f'Cfullonoff^-1: {Cinvfullonoff}')
        print(f'Csimpmax^-1:\n{0.5 * e**2 * Cinvsimpmax / h / self.escale}')
        print(f'Csimpmin^-1:\n{0.5 * e**2 * Cinvsimpmin / h / self.escale}')
        print(f'Csimponoff^-1: {Cinvsimponoff}')

        print(f'Cijmax:\n{Cijmax / self.scale}')
        print(f'Cijmin:\n{Cijmin / self.scale}')
        print(f'Cijonoff: {Cijonoff}')

        # Write the figure to file
        myplt.write_fig_to_file(fname)

    def plot_coupling_cmatrix(self, fname, parasitic=False,
                              xrot=None, ylim=None, vg_start_idx=None):
        """
        Computes and plots the coupling matrix in the full parasitic and
        simplified limits of the coupler
        """
        # Physical constants to convert from mks to MHz
        e = 1.602176634e-19
        h = 6.62607015e-34

        # Get the full C-matrix vs. Vg
        Vg  = np.array([Cval[0] for Cval in self.data.values()])

        # Common plotting functions to both
        myplt = plt.MPLPlotWrapper()
        idx = 0
        mrks = myplt.marker_cycle
        mlen = len(mrks)

        # Full parasitic case
        Cinv = self.get_coupling_matrix(parasitic=parasitic)
        shift = 0 if parasitic else 1

        # Plot the matrix elements vs. Vg
        for i in range(self.Nt - shift):
            for j in range(self.Nt - shift):
                print(f'C[{i}, {j}]')
                if vg_start_idx is not None:
                    myplt.plot(Vg[vg_start_idx:],
                      0.5 * e**2 * Cinv[vg_start_idx:, i, j] / h / self.escale,
                               marker=mrks[idx % mlen], ls='-',
                        label=r'$C^{-1}_{%d%d}$' % (i + 1, j + 1))
                else:
                    myplt.plot(Vg, 0.5 * e**2 * Cinv[:, i, j] / h / self.escale,
                               marker=mrks[idx % mlen], ls='-',
                        label=r'$C^{-1}_{%d%d}$' % (i + 1, j + 1))
                idx += 1

        # Final plotting command common to both
        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$\frac{1}{2}e^2C^{-1}_{ij}$ [MHz]'
        if xrot is not None:
            myplt.set_xaxis_rot(xrot)
        if ylim is not None:
            myplt.ylim = ylim
        myplt.set_leg_outside()

        # Write the figure to file
        myplt.write_fig_to_file(fname)

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
        myplt.plot(Vg, np.abs(Cij) / self.scale, 'o-')
        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$|C_{%d%d}|$ [fF]' % (i + 1, j + 1)

        # Write the figure to file
        myplt.write_fig_to_file(fname)

    def plot_c_all_vs_vg(self, fname, ylim=None, xrot=None, vg_start_idx=None):
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
                if vg_start_idx is not None:
                    myplt.plot(Vg[vg_start_idx:],
                               Cij[vg_start_idx:, i, j] / self.scale,
                               marker=mrks[idx % mlen], ls='-',
                        label=r'$C_{%d%d}$' % (i + 1, j + 1))
                else:
                    myplt.plot(Vg, Cij[:, i, j] / self.scale,
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

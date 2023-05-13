# -*- coding: utf-8 -*-
"""
Tools to analyze the admittance matrix data
"""

import plotting as plt
import numpy as np
import scipy.interpolate
import re as regex
from scipy.optimize import curve_fit

class Admittance(object):
    """
    Class for handling the admittance data
    """
    def __init__(self, *args, **kwargs):

        # Default class conditions
        self.all_data_loaded = False
        self.units = r'$C$ [fF]'
        self.cscale = 1e-15
        self.gscale = 1e6
        self.fscale = 1e9
        self.tscale = 1e-9
        self.symmetrize = False

        # Update the arguments and keyword arguments
        self.__dict__.update(locals())
        for k, v in kwargs.items():
            setattr(self, k, v)


    def __del__(self):
        pass

    def load_all_data_from_file_fourier(self, fname, skip_header=5):
        """
        Loads the admittance data from a single file using the header
        data to extract the gate voltages, times, then apply a Fourier transform
        to the time domain data to recover the frequency content
        """
        # Read the data and order into Vg-labeled dictionary
        data       = np.genfromtxt(fname, dtype=complex,
                                   skip_header=skip_header).T
        # Zero out the nans
        data[np.isnan(data)] = 0
        self.time  = np.unique(data[4].real)
        Ntpts      = self.time.size
        Vg_uniq    = np.unique(data[0].real)[::-1]
        self.Vg    = Vg_uniq
        Nvg        = Vg_uniq.size
        Y          = data[6:]
        self.Nt = int(np.sqrt(Y.shape[0]))
        self.data = {}
        Y = np.reshape(Y, [self.Nt * self.Nt, Nvg, self.Nt, Ntpts])
        for vidx, vg in enumerate(Vg_uniq):
            for tidx, t in enumerate(self.time):
                Y_key = f'y_vg_{vg}_t_{t/self.tscale}ns'.replace('.', 'p')
                Yin = np.zeros(self.Nt * self.Nt, dtype=Y.dtype)
                for i in range(self.Nt):
                    # Yin += Y[vidx, i, fidx, :]
                    # Yin += Y[:, fidx, i, vidx]
                    Yin += Y[:, vidx, i, tidx]
                self.data[Y_key] = [vg, Yin.reshape([self.Nt, self.Nt])]

        self.all_data_loaded = True

    def load_all_data_from_iv_files(self, ifname, vfname, skip_header=5,
                                    debug=False):
        """
        Loads the admittance data from a pair of files using the header
        data to extract the gate voltages, times, then apply a Fourier transform
        to the time domain data to recover the frequency content
        """
        # Read the data and order into Vg-labeled dictionary
        idata       = np.genfromtxt(ifname, dtype=complex,
                                   skip_header=skip_header).T
        vdata       = np.genfromtxt(vfname, skip_header=skip_header).T

        # Zero out the nans
        idata[np.isnan(idata)] = 0
        vdata[np.isnan(vdata)] = 0
        self.freq  = np.unique(idata[4].real)
        Nf         = self.freq.size
        Vg_uniq    = np.unique(idata[0].real)[::-1]
        self.Vg    = Vg_uniq
        Nvg        = Vg_uniq.size
        Ynum       = idata[6:]
        Yden       = vdata[6:]
        self.Nt = Ynum.shape[0]
        print(f'Nt: {self.Nt}')
        self.data = {}
        Ynum = np.reshape(Ynum, [self.Nt, Nvg, self.Nt, Nf])
        Yden = np.reshape(Yden, [self.Nt, Nvg, self.Nt, Nf])
        for vidx, vg in enumerate(Vg_uniq):
            for fidx, f in enumerate(self.freq):
                Y_key = f'y_vg_{vg}_f_{f/self.fscale}GHz'.replace('.', 'p')
                Yin_num = np.zeros(self.Nt * self.Nt, dtype=Ynum.dtype)
                Yin_den = np.zeros(self.Nt * self.Nt, dtype=Yden.dtype)
                Yin = np.zeros([self.Nt, self.Nt], dtype=Ynum.dtype)
                for i in range(self.Nt):
                    for j in range(self.Nt):
                        # Yin_num = Ynum[i, vidx, j, fidx]
                        # for k in range(self.Nt):
                        #     Yin_den = Yden[i, vidx, k, fidx]
                        for k in range(self.Nt):
                            Yin_num = Ynum[i, vidx, k, fidx]
                            Yin_den = Yden[k, vidx, j, fidx]
                            if not np.isclose(Yin_den, 0.):
                                if debug:
                                    s1 = f'I / V[{i}, {j}]'
                                    s2 = f'{Yin_num:.3e} / {Yin_den:.1e}'
                                    print(f'{s1}: {s2}')
                                Yin[i, j] = Yin_num / Yin_den
                self.data[Y_key] = [vg, Yin]

        self.all_data_loaded = True

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
        Vg_uniq    = np.unique(data[0].real)[::-1]
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

    def yij_to_cij_vs_vg(self, i, j):
        """
        Convert the frequency-dependent Y_ij to C_ij(V_g)
        """
        # Read the data into an array
        Cij = np.zeros(self.Vg.size)
        def fitfun(freq, c0, c1):
            return c0 * freq + c1
        for vidx, Vg in enumerate(self.Vg):
            # Get the data at Vg
            Y          = self.group_by_frequency_at_vg(Vg)
            Yij        = np.array([YY[i, j] for YY in Y])
            # popt, pcov = curve_fit(fitfun, 2*np.pi*self.freq, Yij.imag)
            Cij[vidx]  = np.imag(np.gradient(Yij, 2 * np.pi * self.freq)).mean()

        return Cij

    def yij_to_gij_vs_vg(self, i, j):
        """
        Convert the frequency-dependent Y_ij to G_ij(V_g)
        """
        # Read the data into an array
        Gij = np.zeros(self.Vg.size)
        def fitfun(freq, g0, g1):
            return g0 * freq + g1

        for vidx, Vg in enumerate(self.Vg):
            # Get the data at Vg
            Y          = self.group_by_frequency_at_vg(Vg)
            Yij        = np.array([YY[i, j] for YY in Y])
            popt, pcov = curve_fit(fitfun, self.freq, Yij.real)
            Gij[vidx]  = popt[1] # popt[0] * (self.freq[1] - self.freq[0])

        return Gij

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

    def plot_c_all_vs_vg(self, fname, ylim=None, xrot=None):
        """
        Plot the full capacitance matrix as a function of gate voltage
        """
        # Initialize the figure and axes
        myplt = plt.MPLPlotWrapper()
        idx = 0
        mrks = myplt.marker_cycle
        mlen = len(mrks)
        for i in range(self.Nt):
            for j in range(self.Nt):
                print(f'C[{i}, {j}]')
                Cij = self.yij_to_cij_vs_vg(i, j)
                myplt.plot(self.Vg, Cij / self.cscale,
                           marker=mrks[idx % mlen], ls='-',
                    label=r'$C_{%d%d}$' % (i + 1, j + 1))
                idx += 1

        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$C_{ij}$ [fF]'
        myplt.ylim = ylim
        if xrot is not None:
            myplt.set_xaxis_rot(xrot)
        myplt.set_leg_outside()

        # Write the figure to file
        myplt.write_fig_to_file(fname)


    def plot_g_all_vs_vg(self, fname, ylim=None, xrot=None):
        """
        Plot the full conductance matrix as a function of gate voltage
        """
        # Initialize the figure and axes
        scale_dict = {'G' : 1e9, 'M' : 1e6, 'k' : 1e3, '' : 1.,
                'm' : 1e-3, r'$\mu$' : 1e-6}
        scale_str = \
        list(scale_dict.keys())[list(scale_dict.values()).index(self.gscale)]
        myplt = plt.MPLPlotWrapper()
        idx = 0
        mrks = myplt.marker_cycle
        mlen = len(mrks)
        for i in range(self.Nt):
            for j in range(self.Nt):
                print(f'G[{i}, {j}]')
                Gij = self.yij_to_gij_vs_vg(i, j)
                myplt.plot(self.Vg, Gij / self.gscale,
                           marker=mrks[idx % mlen], ls='-',
                    label=r'$G_{%d%d}$' % (i + 1, j + 1))
                idx += 1

        myplt.xlabel = r'$V_g$ [V]'
        myplt.ylabel = r'$G_{ij}$ [%sS]' % scale_str
        myplt.ylim = ylim
        if xrot is not None:
            myplt.set_xaxis_rot(xrot)
        myplt.set_leg_outside()

        # Write the figure to file
        myplt.write_fig_to_file(fname)

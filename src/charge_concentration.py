# -*- coding: utf-8 -*-
"""
Tools to analyze the electron charge concentrations and extract the effect
capacitance widths
"""

import plotting as plt
import numpy as np
import scipy.interpolate
import re as regex

class ChargeConcentration(object):
    """
    Class for handling the charge concentration data
    """
    def __init__(self, *args, **kwargs):

        # Default class conditions
        self.all_data_loaded = False
        self.units = r'log$_{10}$ n[cm$^{-3}]$'
        self.scale = 6

        # Update the arguments and keyword arguments
        self.__dict__.update(locals())
        for k, v in kwargs.items():
            setattr(self, k, v)


    def __del__(self):
        pass

    def load_data_from_file(self, fnames, labels, skip_header=8):
        """
        Loads charge concentration data from a text file exported from a COMSOL
        plot of the log10N data
        """
        # Load the data for each filename
        self.data = {}
        for fname, label in zip(fnames, labels):
            print(f'Loading ({label})-data from {fname} ...')
            self.data[label] = np.genfromtxt(fname, delimiter=',',
                    skip_header=skip_header).T

        self.all_data_loaded = False

    def load_all_data_from_file(self, fname, skip_header=9):
        """
        Loads the charge concentration data from a single file using the header
        data to extract the gate voltages
        """
        # Read the data as a string and pull out the header
        with open(fname, 'r') as fid:
            sdata = fid.read()
            ssplit = sdata.split('\n')
            columns = ssplit[8]
            recmp = regex.compile('Vdc=(.*)')
            Vgs = np.float64(recmp.findall(columns.replace(',', '\n')))

        # Read the data and order into Vg-labeled dictionary
        data = np.genfromtxt(fname, delimiter=',', skip_header=skip_header).T
        x = data[0]
        y = data[1]
        z = data[2:] + self.scale
        print(f'z[0].shape: {z[0].shape}')
        self.data = {}
        print(f'Vgs: {Vgs}')
        for vidx, vg in enumerate(Vgs):
            vg_key = f'ng_vg_{vg}_charge'.replace('.', 'p')
            self.data[vg_key] = [x, y, z[vidx]]

        self.all_data_loaded = True

    def plot_concentrations(self, fig_dir='figs/',
                            xyunits={'x' : r'[$\mu$m]', 'y' : '[nm]',
                                'xv' : 1e-6, 'yv' : 1e-9},
                            binary_threshold=None, plot_option='scatter',
                            xlim=None, ylim=None):
        """
        Plots all concentration data using the labels as filenames
        """
        # Two tone color for values < binary_threshold
        if binary_threshold is not None:
            # Iterate over all entries in the data
            idx = 0
            for key, val in self.data.items():
                # Annonate with the key label
                ksplit = key.split('_')
                t1 = ksplit[1].capitalize()
                tt = r'$%s_{%s}$' % (t1[0], t1[1])
                t2 = ksplit[2].replace('p', '.')
                tstr = f'{tt} = {t2} V'

                # Instantiate the plotting class object
                myplt = plt.MPLPlotWrapper()
                x, y, z = val[0], val[1], val[2]

                # Rescale the x and y values
                if self.all_data_loaded:
                    if idx == 0:
                        x /= xyunits['xv']
                        y /= xyunits['yv']
                else:
                    x /= xyunits['xv']
                    y /= xyunits['yv']
                
                # Crop the data to the xlimits
                if xlim:
                    myplt.xlim = xlim
                    xidx = np.where(
                            np.logical_and(x >= xlim[0], x <= xlim[1]))[0]
                    x = x[xidx]
                    y = y[xidx]
                    z = z[xidx]

                z = z < binary_threshold
                myplt.plot_2d_cmap(x, y, z, norm_type='linear',
                        fname=f'./{fig_dir}{key}_thresholded_concentration.pdf',
                        cmap_str='binary', xstr=r'$x$ %s' % xyunits['x'],
                        ystr=r'$y$ %s' % xyunits['y'], plot_option=plot_option,
                        zlim=(0, 1), tstr=tstr)
                idx += 1

        else:
            idx = 0
            for key, val in self.data.items():
                # Annonate with the key label
                ksplit = key.split('_')
                t1 = ksplit[1].capitalize()
                tt = r'$%s_{%s}$' % (t1[0], t1[1])
                t2 = ksplit[2].replace('p', '.')
                tstr = f'{tt} = {t2} V'

                # Instantiate the plotting class object
                myplt = plt.MPLPlotWrapper()
                x, y, z = val[0], val[1], val[2]

                # Rescale the x and y values
                if self.all_data_loaded:
                    if idx == 0:
                        x /= xyunits['xv']
                        y /= xyunits['yv']
                else:
                    x /= xyunits['xv']
                    y /= xyunits['yv']

                # Crop the data to the xlimits
                if xlim:
                    myplt.xlim = xlim
                    xidx = np.where(
                            np.logical_and(x >= xlim[0], x <= xlim[1]))[0]
                    x = x[xidx]
                    y = y[xidx]
                    z = z[xidx]
                if ylim:
                    myplt.ylim = ylim
                    yidx = np.where(
                            np.logical_and(y >= ylim[0], y <= ylim[1]))[0]
                    x = x[yidx]
                    y = y[yidx]
                    z = z[yidx]

                myplt.plot_2d_cmap(x, y, z, norm_type='linear',
                        fname=f'./{fig_dir}{key}_concentration.pdf',
                        # cmap_str='cividis', xstr=r'$x$ %s' % xyunits['x'],
                        cmap_str='viridis', xstr=r'$x$ %s' % xyunits['x'],
                        ystr=r'$y$ %s' % xyunits['y'], plot_option=plot_option,
                        zlim=(-30, 24), cbar_str=self.units, tstr=tstr)

                myplt.close('all')
                idx += 1

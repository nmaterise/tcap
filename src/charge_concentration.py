# -*- coding: utf-8 -*-
"""
Tools to analyze the electron charge concentrations and extract the effect
capacitance widths
"""

import plotting as plt
import numpy as np
import scipy.interpolate

class ChargeConcentration(object):
    """
    Class for handling the charge concentration data
    """
    def __init__(self, *args, **kwargs):
        pass

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

    def plot_concentrations(self, fig_dir='figs/',
                            xyunits={'x' : r'[$\mu$m]', 'y' : '[nm]',
                                'xv' : 1e-6, 'yv' : 1e-9}):
        """
        Plots all concentration data using the labels as filenames
        """
        # Iterate over all entries in the data
        for key, val in self.data.items():
            # Instantiate the plotting class object
            myplt = plt.MPLPlotWrapper()

            # Interpolate, then send to the function
            print(f'Interpolating {key} data ...')
            x, y, z = val[0], val[1], val[2]
            # Rescale the x and y values
            x /= xyunits['xv']
            y /= xyunits['yv']
            myplt.xlim = [10, 10.1]
            myplt.plot_2d_cmap(x, y, z, norm_type='linear',
                    fname=f'./{fig_dir}{key}_concentration.pdf',
                    cmap_str='jet', xstr=r'$x$ %s' % xyunits['x'],
                    ystr=r'$y$ %s' % xyunits['y'], plot_option='scatter',
                    zlim=(1, 19))

# -*- coding: utf-8 -*-
"""
Plot the charge concentrations
"""

import sys
sys.path.append('./src')
import charge_concentration as cc

def test_plot_charges():
    """
    Plots the charge concentrations for a given gate voltage
    """
    dpath = './data'
    vstrs = ['0', 'n3', 'n3p5', 'n4', 'n4p5', 'n5', 'n15']
    labels = [f'n_vg_{v}' for v in vstrs]
    fnames = [f'{dpath}/surface_density_vg{v}.txt' for v in vstrs]

    # Call the charge concentration class
    mycc = cc.ChargeConcentration()
    mycc.load_data_from_file(fnames, labels)
    mycc.plot_concentrations(fig_dir='figs/', plot_option='tricontourf',
                             xlim=(10, 10.1))

def test_plot_charges_binary():
    """
    Plots the charge concentrations for a given gate voltage
    """
    dpath = './data'
    vstrs = ['0', 'n3', 'n3p5', 'n4', 'n4p5', 'n5', 'n15']
    labels = [f'n_vg_{v}' for v in vstrs]
    fnames = [f'{dpath}/surface_density_vg{v}.txt' for v in vstrs]

    # Call the charge concentration class
    mycc = cc.ChargeConcentration()
    mycc.load_data_from_file(fnames, labels)
    mycc.plot_concentrations(fig_dir='figs/', binary_threshold=5)

def test_plot_all_charges():
    """
    Plots the charge concentrations for a given gate voltage
    """
    dpath = './data'
    fname = f'{dpath}/surface_density_50K_vgall.txt'
    fname = f'{dpath}/surface_density_50K_vgall_230618.txt'

    # Call the charge concentration class
    mycc = cc.ChargeConcentration()
    mycc.load_all_data_from_file(fname)
    mycc.plot_concentrations(fig_dir='figs/', plot_option='tricontourf')

def test_plot_all_charges_binary():
    """
    Plots the charge concentrations for a given gate voltage
    """
    dpath = './data'
    fname = f'{dpath}/surface_density_50K_vgall.txt'

    # Call the charge concentration class
    mycc = cc.ChargeConcentration()
    mycc.load_all_data_from_file(fname)
    mycc.plot_concentrations(fig_dir='figs/', plot_option='tricontourf',
            binary_threshold=0)

if __name__ == '__main__':
    # Call the function above to test the source code
    # test_plot_charges()
    # test_plot_charges_binary()
    test_plot_all_charges()
    # test_plot_all_charges_binary()

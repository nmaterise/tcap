# -*- coding: utf-8 -*-
"""
Calls the testing codes with appropriate conditions to produce the figures in
the main text of the paper https://arxiv.org/abs/2212.04598
"""

import sys
sys.path.append('./src')
import datetime
import charge_concentration as cc
import capacitance_matrix as cm
import admittance_matrix as am

dstr = datetime.datetime.today().strftime('%y%m%d')

def test_plot_all_charges():
    """
    Plots the charge concentrations for a given gate voltage
    """
    dpath = './data'
    fname = f'{dpath}/surface_density_50K_vgall_230618.txt'

    # Call the charge concentration class
    mycc = cc.ChargeConcentration()
    mycc.load_all_data_from_file(fname)
    mycc.plot_concentrations(fig_dir='figs/', plot_option='tricontourf')

def test_plot_all_admittance_data():
    """
    Tests the plotting of the Y data vs. Vg
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_freq_vgall_230618.txt'
    myy.load_all_data_from_file(filename)
    myy.cscale = 1e-15
    myy.gscale = 1e-6

    fname = f'./figs/ymatrix_c_all_vs_vg_{dstr}.pdf'
    myy.plot_c_all_vs_vg(fname, xrot=45.) # , ylim=(-15, 0))
    fname = f'./figs/ymatrix_g_all_vs_vg_{dstr}.pdf'
    myy.plot_g_all_vs_vg(fname, xrot=45.) # , ylim=(None, 50))

def test_plot_coupling_cmatrix():
    """
    Tests the calculation of the coupling matrix in the full parasitic and
    simplified versions
    """
    # Instantiate the class and read the data from file
    myc = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_freq_vgall_230618.txt'
    myc.load_all_data_from_file(filename)
    ylim = (-5, 225)
    fname = './figs/simplified_coupling_cmatrix_from_y_vs_vg.pdf'
    myc.plot_coupling_cmatrix(fname, parasitic=False, xrot=45, ylim=ylim)
    fname = './figs/parasitic_coupling_cmatrix_from_y_vs_vg.pdf'
    myc.plot_coupling_cmatrix(fname, parasitic=True, xrot=45, ylim=ylim)

def test_plot_cij_parasitic_simplified_compare():
    """
    Plots a single matrix element of C^-1 from the parasitic and simplified
    cases, i.e. compare the calculation of the exchange rate
    """
    # Setup the class, read the data, plot the matrix elements
    myc = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_freq_vgall_230618.txt'
    myc.load_all_data_from_file(filename)
    fname = './figs/c12_coupling_cmatrix_comparison_from_y_vs_vg.pdf'
    ylim = (-1, 25)
    myc.plot_cij_parasitic_simplified_compare(0, 1, fname, xrot=45, ylim=ylim)

if __name__ == '__main__':
    # 1.) Functions from src/charge_concentration.py -- figure 3
    test_plot_all_charges()

    # 2.) Functions from src/admittance.py
    ## Plot the capacitance and conductance matrices -- figure 4 
    test_plot_all_admittance_data()

    ## Coupling matrix calculations -- figure 7
    test_plot_coupling_cmatrix()
    test_plot_cij_parasitic_simplified_compare()

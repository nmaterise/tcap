# -*- coding: utf-8 -*-
"""
Plot the capacitances vs. gate voltage
"""

import sys
sys.path.append('./src')
import capacitance_matrix as cm

def test_load_capacitance_data():
    """
    Tests the reading of the data for correct ordering, serves as a check before
    implementing terminal-wise averaging
    """
    # Instantiate the class and read the data from file
    myc = cm.MaxwellCapacitance(symmetrize=True)
    filename = './data/cmatrix_50K_vgall.txt'
    myc.load_all_data_from_file(filename)
    C = myc.data
    for k, v in C.items():
        print(f'C[{k}]:\n{v[1]}')

def test_plot_capacitance_data():
    """
    Tests the plotting of the C data vs. Vg
    """
    # Instantiate the class and read the data from file
    myc = cm.MaxwellCapacitance()
    filename = './data/cmatrix_50K_vgall.txt'
    myc.load_all_data_from_file(filename)
    for i in range(myc.Nterm):
        for j in range(i, myc.Nterm):
            fname = f'./figs/cmatrix_c{i+1}{j+1}_vs_vg.pdf'
            myc.plot_cij_vs_vg(i, j, fname)

def test_plot_all_capacitance_data():
    """
    Tests the plotting of the C data vs. Vg
    """
    # Instantiate the class and read the data from file
    myc = cm.MaxwellCapacitance(symmetrize=False)
    filename = './data/cmatrix_50K_vgall.txt'
    myc.load_all_data_from_file(filename)
    fname = './figs/cmatrix_c_all_vs_vg.pdf'
    # myc.plot_c_all_vs_vg(fname, ylim=[-160, None], xrot=45)
    myc.plot_c_all_vs_vg(fname, xrot=45, vg_start_idx=1)

def test_plot_coupling_cmatrix():
    """
    Tests the calculation of the coupling matrix in the full parasitic and
    simplified versions
    """
    # Instantiate the class and read the data from file
    myc = cm.MaxwellCapacitance(symmetrize=False)
    filename = './data/cmatrix_50K_vgall.txt'
    myc.load_all_data_from_file(filename)
    ylim = [5, 200]
    fname = './figs/simplified_coupling_cmatrix_vs_vg.pdf'
    myc.plot_coupling_cmatrix(fname, parasitic=False, xrot=45, ylim=ylim,
            vg_start_idx=1)
    fname = './figs/parasitic_coupling_cmatrix_vs_vg.pdf'
    myc.plot_coupling_cmatrix(fname, parasitic=True, xrot=45, ylim=ylim,
            vg_start_idx=1)

def test_plot_cij_parasitic_simplified_compare():
    """
    Plots a single matrix element of C^-1 from the parasitic and simplified
    cases, i.e. compare the calculation of the exchange rate
    """
    # Setup the class, read the data, plot the matrix elements
    myc = cm.MaxwellCapacitance(symmetrize=False)
    filename = './data/cmatrix_50K_vgall.txt'
    myc.load_all_data_from_file(filename)
    fname = './figs/c12_coupling_cmatrix_comparison_vs_vg.pdf'
    myc.plot_cij_parasitic_simplified_compare(0, 1, fname, xrot=45, ylim=None,
            vg_start_idx=1)

if __name__ == '__main__':
    # Call the function above to test the source code
    # test_load_capacitance_data()
    # test_plot_capacitance_data()
    # test_plot_all_capacitance_data()
    # test_plot_coupling_cmatrix()
    test_plot_cij_parasitic_simplified_compare()

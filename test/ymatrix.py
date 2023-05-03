# -*- coding: utf-8 -*-
"""
Plot the admittances vs. gate voltage
"""

import sys
sys.path.append('./src')
import admittance_matrix as am

def test_load_admittance_data():
    """
    Tests the reading of the data for correct ordering, serves as a check before
    implementing terminal-wise averaging
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_vgall.txt'
    myy.load_all_data_from_file(filename)
    Y = myy.data
    for k, v in Y.items():
        print(f'Y[{k}]:\n{v[1]}')

def test_plot_admittance_data():
    """
    Tests the plotting of the Y data vs. Vg
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_vgall.txt'
    myy.load_all_data_from_file(filename)
    Vg = -4.5
    Vgstr = f'{Vg}'.replace('.', 'p')
    for i in range(myy.Nt):
        for j in range(myy.Nt):
            fname = f'./figs/ymatrix_y{i+1}{j+1}_vg_{Vgstr}V.pdf'
            myy.plot_yij_vs_f_at_vg(i, j, Vg, fname)

def test_plot_admittance_data_vs_vg():
    """
    Tests the plotting of the Y data vs. Vg
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_vgall.txt'
    myy.load_all_data_from_file(filename)
    real_imag = 'imag'
    for i in range(myy.Nt):
        for j in range(myy.Nt):
            fname = f'./figs/ymatrix_y{real_imag[0:2]}{i+1}{j+1}_vgall.pdf'
            myy.plot_yij_vs_f_vs_vg(i, j, fname, real_imag=real_imag)

def test_plot_all_admittance_data():
    """
    Tests the plotting of the Y data vs. Vg
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=True)
    filename = './data/ymatrix_50K_vgall.txt'
    myy.load_all_data_from_file(filename)
    fname = './figs/ymatrix_c_all_vs_vg.pdf'
    myy.plot_c_all_vs_vg(fname)

if __name__ == '__main__':
    # Call the function above to test the source code
    # test_load_admittance_data()
    # test_plot_admittance_data()
    # test_plot_admittance_data()
    test_plot_admittance_data_vs_vg()

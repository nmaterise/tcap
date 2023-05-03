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
    myy = am.Admittance(symmetrize=True)
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
    myy = am.Admittance()
    filename = './data/ymatrix_50K_vgall.txt'
    myy.load_all_data_from_file(filename)
    for i in range(myy.Nterm):
        for j in range(i, myy.Nterm):
            fname = f'./figs/ymatrix_c{i+1}{j+1}_vs_vg.pdf'
            myy.plot_cij_vs_vg(i, j, fname)

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
    test_load_admittance_data()
    # test_plot_admittance_data()
    # test_plot_all_admittance_data()

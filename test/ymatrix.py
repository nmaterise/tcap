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

def test_load_admittance_data_iv():
    """
    Tests reading the data from separate files for I, V; with V being the
    differential voltage and I the terminal current from a harmonic perturbation
    calculation
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=False)
    ifname = './data/ivector_50K_vgall.txt'
    vfname = './data/vvector_50K_vgall.txt'
    myy.load_all_data_from_iv_files(ifname, vfname, debug=True)
    Y = myy.data
    # for k, v in Y.items():
    #     print(f'Y[{k}]:\n{v[1]}')

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
    filename = './data/ymatrix_50K_freq_eq4_vgall.txt'
    # filename = './data/ymatrix_50K_stat_freq_vgall.txt'
    myy.load_all_data_from_file(filename)
    real_imag = 'imag'
    for i in range(myy.Nt):
        for j in range(myy.Nt):
            fname = f'./figs/ymatrix_y{real_imag[0:2]}{i+1}{j+1}_vgall.pdf'
            myy.plot_yij_vs_f_vs_vg(i, j, fname, real_imag=real_imag)
    real_imag = 'real'
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
    filename = './data/ymatrix_50K_stat_freq_vgall.txt'
    filename = './data/ymatrix_50K_stationary_vgall.txt'
    filename = './data/ymatrix_50K_vgall.txt'
    # filename = './data/ymatrix_50K_freq_eq_vgall.txt'
    filename = './data/ymatrix_50K_freq_eq2_vgall.txt'
    filename = './data/ymatrix_50K_freq_eq3_vgall.txt'
    filename = './data/ymatrix_50K_freq_eq4_vgall.txt'
    myy.load_all_data_from_file(filename)
    fname = './figs/ymatrix_c_all_vs_vg.pdf'
    myy.cscale = 1e-18
    myy.plot_c_all_vs_vg(fname, ylim=[-125, 100], xrot=45.)

def test_plot_all_admittance_data_iv():
    """
    Tests the plotting of the Y data vs. Vg
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=True)
    filename = './data/ymatrix_50K_vgall.txt'
    filename = './data/ymatrix_50K_stat_freq_vgall.txt'
    filename = './data/ymatrix_50K_stationary_vgall.txt'
    filename = './data/ymatrix_50K_vgall.txt'
    # filename = './data/ymatrix_50K_freq_eq_vgall.txt'
    # filename = './data/ymatrix_50K_freq_eq2_vgall.txt'
    ifname = './data/ivector_50K_vgall.txt'
    vfname = './data/vvector_50K_vgall.txt'
    myy.load_all_data_from_iv_files(ifname, vfname)

    # fname = './figs/ymatrix_iv_c_all_vs_vg.pdf'
    # myy.plot_c_all_vs_vg(fname, ylim=[-125, 100], xrot=45.)
    fname = './figs/ymatrix_iv_g_all_vs_vg.pdf'
    myy.gscale = 1e-6
    myy.plot_g_all_vs_vg(fname, xrot=45)

def test_plot_admittance_data_vs_vg_iv():
    """
    Tests the plotting of the Y data vs. Vg
    """
    # Instantiate the class and read the data from file
    myy = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_vgall.txt'
    # filename = './data/ymatrix_50K_stat_freq_vgall.txt'
    ifname = './data/ivector_50K_vgall.txt'
    vfname = './data/vvector_50K_vgall.txt'
    myy.load_all_data_from_iv_files(ifname, vfname)
    real_imag = 'imag'
    for i in range(myy.Nt):
        for j in range(myy.Nt):
            fname = f'./figs/ymatrix_y{real_imag[0:2]}{i+1}{j+1}_vgall.pdf'
            myy.plot_yij_vs_f_vs_vg(i, j, fname, real_imag=real_imag)
    # real_imag = 'real'
    # for i in range(myy.Nt):
    #     for j in range(myy.Nt):
    #         fname = f'./figs/ymatrix_y{real_imag[0:2]}{i+1}{j+1}_vgall.pdf'
    #         myy.plot_yij_vs_f_vs_vg(i, j, fname, real_imag=real_imag)

if __name__ == '__main__':
    # Call the function above to test the source code
    ## Y_ij computed with I / V in COMSOL
    # test_load_admittance_data()
    # test_plot_admittance_data()
    # test_plot_admittance_data_vs_vg()
    # test_plot_admittance_data()
    # test_plot_all_admittance_data()

    ## I, V computed in COMSOL and ratios calculated with Python
    # test_load_admittance_data_iv()
    test_plot_all_admittance_data_iv()
    # test_plot_admittance_data_vs_vg_iv()

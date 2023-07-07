# -*- coding: utf-8 -*-
"""
Plot the admittances vs. gate voltage
"""

import sys
import datetime
sys.path.append('./src')
import admittance_matrix as am
dstr = datetime.datetime.today().strftime('%y%m%d')

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
    myy = am.Admittance(symmetrize=False)
    filename = './data/ymatrix_50K_vgall.txt'
    filename = './data/ymatrix_50K_stat_freq_vgall.txt'
    filename = './data/ymatrix_50K_stationary_vgall.txt'
    filename = './data/ymatrix_50K_vgall.txt'
    # filename = './data/ymatrix_50K_freq_eq_vgall.txt'
    filename = './data/ymatrix_50K_freq_eq2_vgall.txt'
    filename = './data/ymatrix_50K_freq_eq3_vgall.txt'
    filename = './data/ymatrix_50K_freq_eq4_vgall.txt'
    filename = './data/ymatrix_50K_freq_vgall.txt'
    filename = './data/ymatrix_50K_freq_vgall_230618.txt'
    myy.load_all_data_from_file(filename)
    myy.cscale = 1e-15
    myy.gscale = 1e-6

    fname = f'./figs/ymatrix_c_all_vs_vg_{dstr}.pdf'
    myy.plot_c_all_vs_vg(fname, xrot=45.) # , ylim=(-15, 0))
    fname = f'./figs/ymatrix_g_all_vs_vg_{dstr}.pdf'
    myy.plot_g_all_vs_vg(fname, xrot=45.) # , ylim=(None, 50))

    fname = f'./figs/ymatrix_c_all_vs_vg_inset_{dstr}.pdf'
    myy.plot_c_all_vs_vg(fname, xrot=45., ylim=(-15, 0))
    fname = f'./figs/ymatrix_g_all_vs_vg_inset_{dstr}.pdf'
    myy.plot_g_all_vs_vg(fname, xrot=45., ylim=(-0.75, 0.75))

    # for i in range(myy.Nt):
    #     for j in range(myy.Nt):
    #         fname = f'./figs/g_{i+1}{j+1}_{dstr}.pdf'
    #         myy.plot_yij_vs_f_vs_vg(i, j, fname, real_imag='real')
    # for i in range(myy.Nt):
    #     for j in range(myy.Nt):
    #         fname = f'./figs/c_{i+1}{j+1}_{dstr}.pdf'
    #         myy.plot_yij_vs_f_vs_vg(i, j, fname, real_imag='imag')

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

    fname = './figs/ymatrix_iv_c_all_vs_vg_inset.pdf'
    myy.cscale = 1e-15
    myy.plot_c_all_vs_vg(fname, xrot=45., ylim=(0, -20))
    fname = './figs/ymatrix_iv_g_all_vs_vg_inset.pdf'
    myy.gscale = 1e-15
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
    # Call the function above to test the source code
    ## Y_ij computed with I / V in COMSOL
    test_plot_all_admittance_data()

    # Coupling matrix calculations
    test_plot_coupling_cmatrix()
    test_plot_cij_parasitic_simplified_compare()

# -*- coding: utf-8 -*-
"""
Simplified Csd calculation from charges on 2DEG edges

"""

import numpy as np
import post_proc_tools as ppt
from scipy.optimize import curve_fit
from matplotlib import cm

def load_qsd_vsd_data(fname, vg):
    """
    Loads data from the Qsd vs. Vsd COMSOL plot export
    """

    # Read the data from file
    qdata = np.genfromtxt(fname, skip_header=8).T

    # Extract the vsd values
    vsd = np.unique(qdata[0])

    # Extract the charges
    qsd = qdata[1].reshape([vsd.size, vg.size])

    return vsd, qsd


def load_n_data(fname):
    """
    Loads the electron concentration data from txt file
    """

    # Read the data from the txt file
    ndata = np.genfromtxt(fname, skip_header=8).T

    # Find the voltages
    vg = np.flip(np.unique(ndata[0]))

    # Extract the electron concentration
    n = ndata[1, 0:vg.size]
    
    return n


def fit_csd_from_vsd_qsd(vsd, vg, qsd):
    """
    Computes the fits for the Csd capacitance from the Qsd vs. Vsd data
    """

    # Allocate memory for the fits and variances
    Csd    = np.zeros(vg.size)
    Csderr = np.zeros(vg.size)

    # Linear fit function
    def fitfun(x, a, b):
        return a * x + b

    # Iterate over all vg points
    for i in range(vg.size):
    
        ## Fit the data to a line
        copts, ccovs = curve_fit(fitfun, vsd, qsd[:,i])

        ## Extract the capacitances
        Csd[i] = copts[0]

        ## Extract the standard deviations
        Csderr[i] = np.sqrt(np.diag(ccovs))[0]


    return Csd, Csderr


def plot_csd_vg(fname, vg, Csd, Csderr, n=None):
    """
    Plots the Csd values for a given pair of Vg, Csd
    """

    # Setup the figure
    fsize = 20
    fig, ax = ppt.init_subplots(fsize=fsize, tight_layout=True)

    # Plot the data
    C0 = 1e15
    ax.plot(vg, C0 * Csd, '.-', color='tab:blue', label=r'$C_{sd}$')
    # ppt.stderr_fill(ax, vg, C0 * Csd, C0 * Csderr, 'tab:blue')
    ax.set_xlabel(r'$V_g$ [V]', fontsize=fsize)
    ax.set_ylabel(r'$C_{sd}$ [fF]', fontsize=fsize)
    ax.set_yscale('log')
    ppt.set_xaxis_rot(ax, 45)
    #ax.set_ylim([1e-3, 1e-1])

    # Twinx the electron concentration with the capacitance
    if n is not None:
        ax2 = ax.twinx()
        ppt.set_axes_fonts(ax2, fsize)
        ax2.plot(vg, n, '.-', color='tab:red', label=r'$n^{2/3}$')
        ax2.set_ylabel(r'$n^{2/3}$ [cm$^{-2}$]', fontsize=fsize)
        ax2.set_yscale('log')
        
        ppt.set_leg_hdls_lbs(ax, fsize, loc='upper left')
        ppt.set_leg_hdls_lbs(ax2, fsize, loc='center left')
    
    # Save to file
    ppt.write_fig_to_file(fig, fname, fsize=fsize)
    

def plot_csd_vg_wrapper(fname, vg=[0, -1, -10], nfname=None):
    """
    Wrapper on the above functions
    """

    # Convert the filename to a figure filename
    fstrip = fname.split('/')[-1].split('.')[0]
    figname = 'figs/%s_csd_vg.pdf' % fstrip

    # Convert vg to numpy array
    vg = np.asarray(vg)

    # Read the data and fit the data
    print('Reading data from file (%s) ...' % fname)
    vsd, qsd = load_qsd_vsd_data(fname, vg)

    print('Fitting the capacitances ...')
    Csd, Csderr = fit_csd_from_vsd_qsd(vsd, vg, qsd)

    # Read the electron concentration data from file
    n = load_n_data(nfname) if nfname else None
    
    # Plot the results
    print('Writing figure to file ...')
    plot_csd_vg(figname, vg, Csd, Csderr, n=n)


def plot_cmatrix(fname, Nterm=3, C0=1e-15, onoff='off',
                 C1=96.6e-15, C2=96.6e-15, escale=1e6):
    """
    Plots the capacitance matrix as a 2D image
    """
    # Read in the data from file
    # C = np.genfromtxt(fname, skip_header=5).reshape([Nterm, Nterm])

    # Hard-code C-matrix for single gate on / off limits
    Con = np.array([[45.7, -32.01, -13.5],
                  [-32.01, 45.7, -13.5],
                  [-13.5, -13.5, 27.0]])

    Coff = np.array([[15.15, -1.42, -13.5],
                  [-1.42, 15.15, -13.5],
                  [-13.5, -13.5, 27.0]])

    Coff  = np.array([[3.260310997364569E-14,
                      -3.2120996782048274E-14,
                      -2.549322730517135E-16],
                      [-3.212099678204827E-14,
                      3.2603109973808003E-14,
                      -2.5493227305168575E-16],
                      [-2.5493227305171346E-16,
                       -2.549322730516857E-16,
                       5.100357874819216E-16]])

    Con = np.array([
    [3.237485719690459E-14,-3.171777236445142E-14,  -4.305302716947267E-16],
    [-3.171777236445146E-14, 3.23748571961925E-14,  -4.3053027169688553E-16],
    [-4.305302716947266E-16,-4.3053027169688513E-16, 8.624851553778127E-16]])

    Coff = np.array([
    [ 1.777312380083867E-15, -1.1541663721693618E-15,-3.968694896706115E-16],
    [-1.1541663721693614E-15, 1.777312380109116E-15, -3.968694896705159E-16],
    [-3.9686948967061135E-16, -3.968694896705159E-16, 7.943995935402962E-16]])


    Con = np.array([
    [3.2280538403058546E-14, -3.1641009953701286E-14, -4.1339875523296384E-16],
    [-3.1641009953701273E-14, 3.2280538404354645E-14, -4.133987552264484E-16],
    [-4.133987552329638E-16, -4.1339875522644873E-16, 8.290742786005986E-16] ])

    Coff = np.array([
    [1.6570721976215746E-15,-1.0935814742645778E-15, -3.373906757388025E-16],
    [-1.0935814742645784E-15, 1.6570721976206136E-15, -3.373906757387148E-16],
    [-3.373906757388023E-16, -3.3739067573871467E-16,  6.757954356913526E-16] ])

    onoff_ratio = Con[0, 1] / Coff[0, 1]
    if onoff == 'on':
        C = np.array([[82.6, -31.7, -50.9],
                      [-31.7, 82.6, -50.9],
                      [-50.9, -50.9, 102.]])
        C = np.copy(Con)
    elif onoff == 'off':
        C = np.array([[52.0, -1.11, -50.9],
                      [-1.11, 52.0, -50.9],
                      [-50.9, -50.9, 102.]])
        C = np.copy(Coff)
    else:
        raise ValueError(f'{onoff} option not recognized.')

    # Compute the parasitic capacitance matrix
    C11 = C1 + abs(C[0, 0]) + abs(C[0, 1]) + abs(C[0, 2])
    C12 = -abs(C[0, 1])
    C13 = -abs(C[0, 2])
    C22 = C2 + abs(C[1, 1]) + abs(C[0, 1]) + abs(C[1, 2])
    C23 = -abs(C[1, 2])
    C33 = C[2, 2] + abs(C[0, 2]) + abs(C[1, 2]) 
    Cfull = np.array([[C11, C12, C13],
                      [C12, C22, C23],
                      [C13, C23, C33]])

    # Compute the inverse capacitance matrix and scale by e^2
    Cs = np.array([[C1 + abs(C[0, 1]), -abs(C[0, 1])],
                [-abs(C[0, 1]), C2 + abs(C[0, 1])]])
    Csinv = np.linalg.inv(Cs)
    Cfullinv = np.linalg.inv(Cfull)
    e = 1.602176634e-19
    h = 6.62607015e-34

    # Convert the matrix elements in MHz
    # CsinvMHz = 0.5* e**2*Csinv*1e15/h/escale
    # CfullinvMHz = 0.5* e**2*Cfullinv*1e15/h/escale
    CsinvMHz = 0.5* e**2*Csinv/h/escale
    CfullinvMHz = 0.5* e**2*Cfullinv/h/escale
 
    np.set_printoptions(precision=2)
    print(f'Ccoupler:\n{C} fF')
    print(f'Csimple:\n{Cs} fF')
    print(f'Csimple^-1:\n{CsinvMHz} MHz')
    print(f'Cfull:\n{Cfull} fF')
    print(f'Cfull^-1:\n{CfullinvMHz} MHz') 
    

    # Call the ppt plotting routine
    fstrip = fname.split('/')[-1].split('.')[0]
    fname_out_cbar = 'figs/%scmatrix_%s_%dx%d_cbar.pdf' \
                     % (fstrip, onoff, Nterm, Nterm)
    fname_out_cbar_full = 'figs/%scinvmatrix_%s_%dx%d_cbar_full.pdf' \
                     % (fstrip, onoff, Nterm, Nterm)
    fname_out_cbar_simp = 'figs/%scinvmatrix_%s_%dx%d_cbar_simp.pdf' \
                     % (fstrip, onoff, 2, 2)

    # Define the limits and plot the results
    x = np.linspace(1, Nterm, Nterm)
    y = x
    cmap = cm.viridis
    # ppt.plot_2d_cmap(x, y, C/C0, fname_out, cbar_str='C [fF]',
    #                  use_imshow=True, show_cbar=False, cmap=cmap,
    #                  zlims=[1e-2, 1.5e2])
    ppt.plot_2d_cmap(x, y, C/C0, fname_out_cbar,
                     cbar_str=r'$C_{%s}$ [fF]' % onoff,
                     use_imshow=True, show_cbar=True, cmap=cmap,
                     # zlims=[1e-2, 1.5e2], nlevels=300)
                     zlims=None, nlevels=300)
    ppt.plot_2d_cmap(x, y, CfullinvMHz, fname_out_cbar_full,
                     cbar_str=r'$C^{-1}$ [MHz]',
                     use_imshow=True, show_cbar=True, cmap=cmap,
                     # zlims=[1e-2, 1.5e2], nlevels=300)
                     zlims=None, nlevels=300)
    x2 = np.linspace(1, 2, 2)
    y2 = x2
    ppt.plot_2d_cmap(x2, y2, CsinvMHz, fname_out_cbar_simp,
                     cbar_str=r'$C^{-1}$ [MHz]',
                     use_imshow=True, show_cbar=True, cmap=cmap,
                     # zlims=[1e-2, 1.5e2], nlevels=300)
                     zlims=None, nlevels=300)

    print(f'on/off ratio (simple): {onoff_ratio}')


if __name__ == '__main__':

    # Plot the capacitance from COMSOL
    fname = \
 'data/qsd_vsd_2deg_1000nm_equil_11vg_100nm_ox_10nm_airgap02_150K_3x3_50nm_nplus_100nm_true_gate_three_gates.txt'
    nfname = \
 'data/n_1000nm_2deg02_ox_10nm_500nm_gh_100nm_gw_11vg_50nm_nplus_three_gate.txt'
    vg = np.hstack((0, -np.linspace(1, 4, 13))) 
    vg = np.hstack((0, -np.linspace(1, 4, 4))) 
    vg = np.hstack((0, -np.linspace(1, 10, 10))) 
    vg = np.hstack((0, -np.linspace(1, 11, 6))) 
    
    # plot_csd_vg_wrapper(fname, vg, nfname=nfname)
    # plot_csd_vg_wrapper(fname, vg, nfname=None)
    # fname1 = 'data/cmatrix_electrostatic_one_three_gate_100nm_500nm_2deg_5nm_air.txt'
    # fname2 = 'data/cmatrix_electrostatic_three_gate_50nm_500nm_2deg_5nm_air.txt'
    # fname3 = 'data/cmatrix_electrostatic_three_gate_100nm_500nm_copper_5nm_air.txt'
    # for fname in [fname1, fname2, fname3]:

    plot_cmatrix('', Nterm=3, C0=1e-15, onoff='off')
    plot_cmatrix('', Nterm=3, C0=1e-15, onoff='on')

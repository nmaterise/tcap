# -*- coding: utf-8 -*-
"""
Simplified Csd calculation from charges on 2DEG edges

"""

import numpy as np
import post_proc_tools as ppt
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import cm
import re as regex


def load_y_data(fname, Nterm=3):
    """
    Read the admittance matrix data from file
    """
    # Read the text data to start
    rdata = np.genfromtxt(fname, skip_header=5, delimiter='\t', dtype=str)
    rdata = [np.complex128(r.split()) for r in rdata]
    y = []
    freqs = []
    for r in rdata:
        if len(r) > Nterm:
            y.append(r[1:])
            freqs.append(r[0])
        else:
            y.append(r)
    freqs = np.asarray(freqs).real
    y = np.asarray(y)
    y = y.reshape([y.shape[0]//Nterm, Nterm, Nterm])

    # Compute Z as the matrix inverse of Y
    z = np.array([np.linalg.inv(yy) for yy in y])

    return freqs, y, z


def plot_yij_vs_f(f, yij, fname_out, fsdict={'scale' : 1e9, 'units' : 'GHz'}):
    """
    Plots the real and imaginary components of the admittance matrix as a
    function of the frequency f
    """
    # Setup the figure and run with it from there    
    fsize = 20; lsize = 14
    fscale = fsdict['scale']
    funits = fsdict['units']
    fig, ax = ppt.init_subplots(fsize=fsize)
    Cij = np.mean(np.imag(np.gradient(yij, 2*np.pi*f)))
    ax.plot(f/fscale, yij.real, '.-', label=r'Re $Y$')
    ax.plot(f/fscale, yij.imag, '.-', label=r'Im $Y$, C=%.2g' % Cij)
    ax.set_xlabel(f'Frequency [{funits}]', fontsize=fsize)
    ax.set_ylabel(f'Admittance Matrix Element [S]', fontsize=fsize)
    ppt.set_axes_fonts(ax, fsize)
    leg = ppt.set_leg_outside(ax, lsize)
    ppt.write_fig_to_file(fig, fname_out, leg=leg, is_leg_outside=True)
    plt.close('all')


def plot_rmatrix(fname, Nterm=3, C0=1e-15, onoff='off',
                 C1=96.6e-15, C2=96.6e-15, escale=1e6):
    """
    Plots the capacitance matrix as a 2D image
    """
    # Read in the data from file
    print(f'fname: {fname}')
    freqs, y, z = load_y_data(fname) 
    np.set_printoptions(precision=2)
    C = np.zeros([Nterm, Nterm])
    G = np.zeros([Nterm, Nterm])
    
    # Extract the C-matrix
    for i in range(Nterm):
        for j in range(Nterm):
            # fname_out = './figs/' + fname.split('.')[0] + f'_{i+1}{j+1}.pdf'
            # plot_yij_vs_f(freqs, y[:,i,j], fname_out)
            C[i, j] = np.mean(np.imag(np.gradient(y[:,i,j],
                            2*np.pi*freqs)))
            G[i, j] = np.mean(y[:, i, j].real)

    R = np.linalg.inv(G)

    print(f'C:\n{C/C0} fF')
    print(f'G:\n{G} S')
    print(f'R:\n{R} Ohms')

    # Compute the parasitic capacitance matrix
    C11 = C1 + C[0, 0] + abs(C[0, 1]) + abs(C[0, 2])
    C12 = -abs(C[0, 1]) 
    C13 = -abs(C[0, 2])
    C22 = C2 + C[1, 1] + abs(C[0, 1]) + abs(C[1, 2])
    C23 = -abs(C[1, 2])
    C33 = C[2, 2] + abs(C[0, 2]) + abs(C[1, 2])
    print(f'C33: {C[2, 2]} + {C[0, 2]} + {C[1, 2]} = {C33}')
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
    CsinvMHz = 0.5* e**2*Csinv/h/escale
    CfullinvMHz = 0.5* e**2*Cfullinv/h/escale

    print(f'Cfull:\n{Cfull/C0} fF')
    print(f'Cs:\n{Cs/C0} fF')
    print(f'Cfullinv:\n{CfullinvMHz} MHz')
    print(f'Csinv:\n{CsinvMHz} MHz')

    # Plot the capacitance and conductance matrices
    x = np.linspace(1, Nterm, Nterm)
    y = x
    # Call the ppt plotting routine
    fstrip = fname.split('/')[-1].split('.')[0]
    fname_out_cbar = 'figs/%scmatrix_%s_%dx%d_cbar.pdf' \
                     % (fstrip, onoff, Nterm, Nterm)
    fname_out_gbar = 'figs/%sgmatrix_%s_%dx%d_cbar.pdf' \
                     % (fstrip, onoff, Nterm, Nterm)
    fname_out_cbar_full = 'figs/%scinvmatrix_%s_%dx%d_cbar_full.pdf' \
                     % (fstrip, onoff, Nterm, Nterm)
    fname_out_cbar_simp = 'figs/%scinvmatrix_%s_%dx%d_cbar_simp.pdf' \
                     % (fstrip, onoff, 2, 2)

    cmap = cm.viridis
    ppt.plot_2d_cmap(x, y, C/C0, fname_out_cbar,
                     cbar_str=r'$C_{\mathrm{%s}}$ [fF]' % onoff,
                     use_imshow=True, show_cbar=True, cmap=cmap,
                     # zlims=[1e-2, 1.5e2], nlevels=300)
                     zlims=None, nlevels=300)
    ppt.plot_2d_cmap(x, y, G, fname_out_gbar,
                     cbar_str=r'G [S]',
                     use_imshow=True, show_cbar=True, cmap=cmap,
                     # zlims=[1e-2, 1.5e2], nlevels=300)
                     zlims=[-1e5, 1e5], nlevels=20000, ndisp=2)
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


if __name__ == '__main__':

    # Plot the admittance data, just to inspect it
    fname = './data/ymatrix_split_gnd_plane_one_gate_conducting_220525.txt'
    fname = './data/ymatrix_split_gnd_plane_one_small_gate_conducting_220527.txt'
    fname = './data/ymatrix_50K_depleted_230617.txt'
    plot_rmatrix(fname, Nterm=3, onoff='depleted')
    fname = './data/ymatrix_RT_none_230617.txt'
    plot_rmatrix(fname, Nterm=3, onoff='none')
    fname = './data/ymatrix_50K_conducting_230617.txt'
    plot_rmatrix(fname, Nterm=3, onoff='conducting')
    # fname = './data/ymatrix_split_gnd_plane_one_small_gate_depleted_220527.txt'
    # plot_rmatrix(fname, Nterm=3, onoff='off')
    # plot_rmatrix('', Nterm=3, C0=1, onoff='on')

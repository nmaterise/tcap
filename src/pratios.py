# -*- coding: utf-8 -*-
"""
Read and normalize the participations in the conducting and depleted limits
"""

import numpy as np

def load_pratios_from_file(fname, normalize=True, Nterm=3):
    """
    Read the text file with the participations
    """
    pdata = np.genfromtxt(fname, delimiter='\t', dtype=str, skip_header=5).T
    pdata = np.array([np.float64(p.split()[1:]) for p in pdata])

    if normalize:
        pdata_out = np.zeros(pdata.shape)
        for idx, p in enumerate(pdata):
            pdata_out[idx, :] = p / sum(p) 
        pdata = np.copy(pdata_out)

    return pdata


def pratios_to_T1s(p, tand, w0):
    """
    Computes the T1 limits for each participation ratio
    """
    T1 = 1. / (w0 * p * tand)
    return T1


def test_load_pratios_from_file():
    """
    Test the function above and print the data
    """
    fname = 'pratios_one_gate_depleted_220527.txt'
    pon = load_pratios_from_file(fname)
    fname = 'pratios_one_gate_conducting_220527.txt'
    poff = load_pratios_from_file(fname)

    np.set_printoptions(formatter={'all':lambda x: f'{x:.2E}'})

    # Loss Tangents:
    tand = np.array([4.8e-5, 4.8e-5, 4.8e-5, 4.8e-5, 5e-5, 4.8e-5])
    w0 = 2*np.pi*5e9
    T1on = pratios_to_T1s(pon[-1], tand, w0)
    T1off = pratios_to_T1s(poff[-1], tand, w0)
    T1on_tot = 1. / sum(1. / T1on)
    T1off_tot = 1. / sum(1. / T1off)

    print(f'tand: {tand}')
    print(f'pon:\n{pon[-1]}')
    print(f'poff:\n{poff[-1]}')
    print(f'T1_on: {T1on / 1e-6} us,\nT1on_tot: {T1on_tot/1e-6:.2E} us')
    print(f'T1_off: {T1off / 1e-6} us,\nT1off_tot: {T1off_tot/1e-6:.2E} us')


if __name__ == '__main__':
    # Test loading and normalizing
    test_load_pratios_from_file()

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:21:01 2017

@author: nmaterise
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

sys.path.append(os.path.abspath(os.path.dirname(__file__)).rsplit(os.sep, 1)[0])

import logging
import pyEPR as epr
from pyEPR import ansys as HFSS
epr.logger.setLevel(logging.DEBUG)

# Choose whether to generate results or not
rerun_hfss = True
hfss_completed = False
optimetrics_only = True

# Specify the HFSS project to be analyzed
design_name = 'billy_design_220508'
project_name = 'Tunable_capacitor-qubits-single-value'
pinfo = epr.ProjectInfo(project_path=r'Z:\\tcap_exp\\Tunable-capacitor', \
                            project_name = \
                            # 'Tunable_capacitor-qubits-single-value', \
                            project_name, \
                            # design_name = 'HFSSDesign1')
                            # design_name = design_name)
                            design_name = 'matthieu_design_201012')

# Describe the junctions in the HFSS desgin
pinfo.junctions['j0'] = {
    'rect': 'qubit0_junction',
    'line': 'qubit0_junction_line',
    'Lj_variable': 'qubit0induct',
    'length': 2.0e-6,
}
pinfo.junctions['j1'] = {
    'rect': 'qubit1_junction',
    'line': 'qubit1_junction_line',
    'Lj_variable': 'qubit1induct',
    'length': 2.0e-6,
}
## Artificial JJ for coupler
pinfo.junctions['j2'] = {
    # 'rect': 'JJ01', # Billy's design
    'rect': 'JJ0',
    # 'line': 'jj01_junction_line', # Billy's design
    'line': 'JJ0_line',
    'Lj_variable': 'q01_tuning_inductance',
    'Cj_variable': 'q01_tuning_capacitance',
    'Rj_variable': 'q01_tuning_resistance',
    # 'length': 1.0e-6,
    'length': 100e-9,
}

pinfo.validate_junction_info()

#use file location path:
HFSS_path=os.getcwd()

full_path=HFSS_path+'\\HFSS\\'+project_name+'.aedt'

HFSS_app=HFSS.HfssApp()
HFSS_desktop=HFSS_app.get_app_desktop()

project=HFSS_desktop.open_project(full_path)

if project==None:
    project=HFSS_desktop.new_project()
    project.save(full_path)
    
# project.save(full_path)
project.make_active()
    
solution_type = 'eigenmode'
if design_name in project.get_design_names():
    if overwrite==True:
        project.delete_design(design_name)
        project.save()
        
        # Setup a driven modal design type
        if solution_type == 'driven_modal':
            EM_design=project.new_dm_design(design_name)
        elif solution_type == 'eigenmode':
            EM_design=project.new_em_design(design_name)
        else:
            raise ValueError(f'Solution type ({solution_type}) not recognized.')
    else:
        EM_design=project.get_design(design_name)
        
else:
    if solution_type == 'driven_modal':
        EM_design=project.new_dm_design(design_name)
    elif solution_type == 'eigenmode':
        EM_design=project.new_em_design(design_name)
    else:
        raise ValueError(f'Solution type ({solution_type}) not recognized.')
        
EM_design.make_active()
model=HFSS.HfssModeler(EM_design)

project=HFSS_desktop.open_project(full_path)

if project==None:
    project=HFSS_desktop.new_project()
    project.save(full_path)
    
project.make_active()
    
if design_name in project.get_design_names():
    if overwrite==True:
        project.delete_design(design_name)
        project.save()
        EM_design=project.new_em_design(design_name)
    else:
        EM_design=project.get_design(design_name)
        
else:
    EM_design=project.new_em_design(design_name)
        
EM_design.make_active()

em_setup = EM_design.create_em_setup(name='Test_EM', 
                                   min_freq_ghz=2, 
                                   n_modes=5, 
                                   max_delta_f=1.0, 
                                   min_converged=1, 
                                   converge_on_real=True)

# time.sleep(5)

# Print the project information objects
print(f'Solution type: {pinfo.design.solution_type}')
print(f'\nObject Names:\n{pinfo.get_all_object_names()}')
print(f'\nVariable Names:\n{pinfo.get_all_variables_names()}')
print(f'\nSetup: {pinfo.design.get_setup_names()}')
print(f'\nNumber of eigenmodes: {pinfo.setup.n_modes}')


# Run analysis
if rerun_hfss:
    if not hfss_completed:
        ## Setup the analysis
        if not optimetrics_only:
            print('Rerunning nominal solution setup ...')
            pinfo.setup.analyze()
        
        ## Setup the optimetrics study
        print('Running Optimetrics ...')
        opti_setup = HFSS.Optimetrics(EM_design)
        # setup_name = pinfo.design.optimetrics.get_setup_names()[0]
         
        ## Run the optimetrics study
        filename = 'tunable_coupler_sweep_params.csv'
        # pinfo.design.optimetrics
        opti_setup.create_setup(variable='',
           swp_params=(),
           swp_type='file',
           filename=filename,
           name=f'Sweep_from_file',
           save_fields=True,
           solve_with_copied_mesh_only=False,
           copy_mesh=False) #solve_setup(setup_name)
        opti_setup.solve()
        
    
    ## Interact with hfss
    eprh = epr.DistributedAnalysis(pinfo)
    print(f'pinfo:\n{pinfo}')
    dfile = None
    dfile, variations = eprh.do_EPR_analysis()
    freqs_qs = eprh.get_ansys_frequencies_all()
    print(f'{freqs_qs}')
    if dfile:
        ## Run the quantum analysis of the HFSS results
        print('Setting up quantum analysis ...')
        epra = epr.QuantumAnalysis(eprh.data_filename)
    
        ## Analyze the quantum results
        print('Performing quantum analysis ...')
        results = epra.analyze_all_variations(cos_trunc=8, fock_trunc=7)

    else:
        raise IOError(f'EPR analysis file ({dfile}) does not exist.')

else:
    ## Setup for quantum analysis
    fname = 'depleted_limit_2020-10-13.npz'
    fname = '2020-10-14 12-04-07.npz'
    fname = '2020-10-19 00-24-50.npz'

    # Latest depleted limit data with fF capacitances
    fname = '2022-0608 10-09-32.npz'
    fname = '2022-11-02 17-08-38.npz'
    tmp_path=f'Z:\\tcav\\data\\data-pyEPR\\Tunable_capacitor-qubits-single-value\\billy_design_220508\\{fname}'
    # epra = epr.QuantumAnalysis(eprh.data_filename)
    epra = epr.QuantumAnalysis(tmp_path)
    
    ## Analyze the quantum results
    results = epra.analyze_all_variations(cos_trunc=8, fock_trunc=7)
    
    print(f'results: {results}')
    for res in results:
        print(res)

# if 0:
#     epr_hfss.verbose = True
#     epr_hfss.do_EPR_analysis()
# 
# if 1:  # Hamiltonian analysis
#     today = datetime.today().strftime('%y%m%d')
#     filename = f'Z:\pyEPR\Tunable_capacitor-qubits\HFSSDesign1\HFSSDesign1_{today}.hdf5'  # epr_hfss.data_filename
#     epr = pyEPR_Analysis(filename)
#     results = epr.analyze_all_variations(cos_trunc=8, fock_trunc=7)
#     # epr.plot_Hresults()
#     c = []
#     c01 = []
#     for variation, result in results.items():
#         # print(variation)
#         # print(result)
#         c.append(float(result['hfss_variables']['_q01_tuning_capacitance'][:-2]))
#         c01.append(result['chi_O1'].iloc[0, 1])
# 
# plt.figure()
# plt.plot(c, c01, '+')
# 
# plt.show()

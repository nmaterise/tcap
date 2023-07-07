# Tunable Capacitor For Superconducting Qubits Using an InAs/InGaAs Heterostructure

## Installation with Anaconda
* Install Anaconda [here](https://www.anaconda.com/) with version 3.8 or higher
* Create a new environment `tcap` with the `requirements.txt`
```bash
conda create --name tcap --file requirements.txt
conda activate tcap
```

## Running example Python code
* Run the code to generate the figures in the main text of
  [arXiv:2212.04598](https://arxiv.org/abs/2212.04598)
```bash
python test/generate_all_figures.py
```

## Directories
* src/          Python source code to post-process COMSOL & HFSS data
* test/         Python drivers on the source code produce figures in the paper
* figs/         Figures produced by the test/ code
* data/         Data used in the test/ code from COMSOL & HFSS

## Running COMSOL files
* To compute the charge concentrations and Y matrices, open
  `src/single_small_gate_2deg_inalas_semi_ec.mph` with COMSOL v6.1 or higher
    - The simulation starts by running the first equilibrium study
    - This study sets the initial conditions for the subsequent frequency sweep
    - The Y matrix evaluation group contains the Y vs. frequency data at various
      dc bias points
    - To export the Y matrix data and the charge concentrations, use the Export
      tab at the bottom of the left-hand-side column

## Running the HFSS files
* Run the python script `src/chi_matrix_extract_tcap.py` from `src/` 

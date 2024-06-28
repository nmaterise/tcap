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
* Install [COMSOL](https://www.comsol.com/) version 6.1 or higher with the
  Semiconductor and AC/DC modules to run the following files
* To compute the charge concentrations and Y matrices, open
  `src/single_small_gate_2deg_inalas_semi_ec.mph` with COMSOL
    - The simulation starts by running the first equilibrium study
    - This study sets the initial conditions for the subsequent frequency sweep
    - The Y matrix evaluation group contains the Y vs. frequency data at various
      dc bias points
    - It takes about 7 minutes on my workstation to run this calculation, it may
      take more or less time to run on your system
    - Run Evaluate on the Global Evaluation labeled Admittance Matrix to update
      the data before exporting the Y matrix to file 
    - To export the Y matrix data and the charge concentrations, use the Export
      tab at the bottom of the left-hand-side column
    - After exporting the data, the `i` entries in the file
      `data/ymatrix_vgall*.txt` need to be changed to `j` and `Inf` replaced
      with `nan` to comply with the file parser
* To compute the participations, the file `src/semi_2deg_.mph` performs the
  integrals for the values in Table 2 of the manuscript

## Running the HFSS files
* Run the python script `src/chi_matrix_extract_tcap.py` from `src/` 
```bash
cd src
python chi_matrix_extract_tcap.py
```

## Latest release
[![DOI](https://zenodo.org/badge/663646184.svg)](https://zenodo.org/badge/latestdoi/663646184)

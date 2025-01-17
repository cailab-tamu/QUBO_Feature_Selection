Running qubo.ipynb

0. Pre-Installation Steps
--------------------------
To set up the environment for running the QUBO notebook, follow these steps:

- Create and Activate the Conda Environment:
    conda env create -f environment.yml
    conda activate dwave

- Install Jupyter Kernel (Linux):
    python -m ipykernel install --user --name dwave

- Setup D-Wave API:
    1. Visit https://cloud.dwavesys.com/leap/ and navigate/find Solver API Token.
    2. In your terminal, configure your D-Wave API client:
        dwave config create
    3. When prompted, paste your Solver API Token and hit Enter.

1. Options for Computing QUBO Matrix
------------------------------------
1. Option 1: Generate QUBO Matrix from MATLAB:
    - Run the script qfeatures_driver.m located in ../../../qubo_fs_matlab/efficient_differentiation/.
    - Copy the generated files genes.csv and qubo_matrix.csv to this local directory.

2. Option 2: Use a Pre-Computed QUBO Matrix:
    - Download the pre-computed matrix from the following link:
      https://drive.google.com/drive/folders/1-2og2FAM0_6e3L2_9C7HnOm9sgfrJIRe?usp=sharing

2. Running the QUBO Notebook
----------------------------
- Open qubo.ipynb in Jupyter Notebook.
- Run all cells in the notebook; the selected features will be stored in filt_df_QA.csv.

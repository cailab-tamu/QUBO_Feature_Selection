
# Setting up Qfeatures - QUBO_Feature_Selection

## 0. Pre-Installation Steps (Python)
To set up the environment for running the QUBO notebook, follow these steps:

### 0.1. Create and Activate the Conda Environment
```bash
conda env create -f environment.yml
conda activate dwave
```

### 0.2. Install the `qfeatures` Package
```bash
pip install -e .
```

### 0.3. Install Jupyter Kernel (Linux)
```bash
python -m ipykernel install --user --name dwave
```

### 0.4. Set up D-Wave API
1. Visit [D-Wave Leap](https://cloud.dwavesys.com/leap/) and find your Solver API Token.
2. In your terminal, configure your D-Wave API client:
   ```bash
   dwave config create
   ```
3. When prompted, paste your Solver API Token and hit Enter.

### 0.5. Add MATLAB Path
Add the following directory to your MATLAB path:
```
C:\Users\ssromerogon\Documents\vscode_working_dir\QUBO_Feature_Selection\qfeatures-src-v0.2_matlab
```

## 1. Options for Computing QUBO Matrix
You have three options for generating the QUBO matrix:

### 1.1. Option 1: Generate QUBO Matrix from MATLAB
- Run the script `qfeatures_driver.m` located in `../../../qubo_fs_matlab/efficient_differentiation/`.
- Copy the generated files `genes.csv` and `qubo_matrix.csv` to the current local directory.

### 1.2. Option 2: Use a Pre-Computed QUBO Matrix
- Download the pre-computed matrix from the following link:
  [Google Drive - Pre-Computed QUBO Matrix](https://drive.google.com/drive/folders/1-2og2FAM0_6e3L2_9C7HnOm9sgfrJIRe?usp=sharing)

### 1.3. Option 3: Run the Stand-Alone Python Notebook
- After installing the `qfeatures` package (`pip install -e .`), run the stand-alone notebook `construct_qubo_python.ipynb` in this directory.
- Download the prepared single-cell data from the following link:
  [Google Drive - Single-Cell Data](https://drive.google.com/drive/folders/1-2og2FAM0_6e3L2_9C7HnOm9sgfrJIRe?usp=sharing)
  File: `Data_hESC_EC_day1_5000g_filtered_feature_bc_matrix_h5.h5ad`

## 2. Running the QUBO Notebook
- Open `qubo_solver_from_MATLAB_Q.ipynb` in Jupyter Notebook.
- Run all cells in the notebook.
- The selected features will be saved in `filt_df_QA.csv`.

### 2.1 Running in matlab 
- Open `qfeatures_driver.m` located in `./qubo_fs_matlab/efficient_differentiation/` and execute it (NOTE: install `Quantum Computing` matlab package)
- Cross validation available in `qfeatures_driver_cross_validation.m` in  `./qubo_fs_matlab/efficient_differentiation/`

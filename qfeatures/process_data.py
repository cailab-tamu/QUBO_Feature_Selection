import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

def compute_pearson_residual(matrix, theta=100):
    """
    Compute clipped Pearson residuals for a given sparse matrix, where rows are features and columns are cells.
    Reference: https://doi.org/10.1101/2020.12.01.405886

    Args:
        matrix (scipy.sparse.csr_matrix or numpy.ndarray): 2D sparse or dense array with rows as features and columns as cells.
    
    Returns:
        scipy.sparse.csr_matrix: Clipped Pearson residuals matrix in sparse format.
    """
    # Ensure the input matrix is sparse
    if not issparse(matrix):
        matrix = csr_matrix(matrix)
    
    # Compute the row and column sums
    row_sums = matrix.sum(axis=1)  # Sum along rows (features)
    col_sums = matrix.sum(axis=0)  # Sum along columns (cells)
    total_sum = matrix.sum()       # Total sum of all elements

    # Compute the expected values under independence assumption
    row_sums_dense = np.array(row_sums).flatten()  # Convert to dense
    col_sums_dense = np.array(col_sums).flatten()  # Convert to dense
    expected = (row_sums_dense[:, None] * col_sums_dense[None, :]) / total_sum

    # Ensure expected values are sparse for consistency
    expected_sparse = csr_matrix(expected)

    # Standard deviation for Pearson residuals
    std_dev = np.sqrt(expected_sparse + (expected_sparse.power(2)) / theta)

    # Compute Pearson residuals
    residuals = (matrix - expected_sparse).multiply(std_dev.power(-1))  # Element-wise division

    # Replace NaN and Inf values (caused by division by zero in std_dev) with 0
    residuals.data[np.isnan(residuals.data)] = 0
    residuals.data[np.isinf(residuals.data)] = 0

    # Clip residuals to ±sqrt(n), where n is the number of columns (cells)
    num_cells = matrix.shape[1]
    clip_value = np.sqrt(num_cells)
    residuals.data = np.clip(residuals.data, -clip_value, clip_value)

    return residuals


def mat_transform(X, y):
    # Apply pearson residuals
    #sc.experimental.pp.normalize_pearson_residuals(adata, theta=100)
    #X = np.array(adata.X.T)
    X = compute_pearson_residual(X)
    X = X.toarray()
    y_reshaped = y.reshape(1, -1) 
    print("Data shape : {X.shape}")
    print("Predictor shape : {y_reshaped.shape}")

    Xy = np.vstack((X, y_reshaped)) 
    Xy = csr_matrix(Xy)

    print("Final shape of Xy: {Xy.shape}")


import scanpy as sc

def mat_transform_sc(adata, obs_key):
    # Apply pearson residuals
    X = adata.X.T
    X = compute_pearson_residual(X)
    X = X.toarray()
    y = np.array(adata.obs[obs_key].values)
    y_reshaped = y.reshape(1, -1) 
    print(X.shape)
    print(y_reshaped.shape)

    Xy = np.vstack((X, y_reshaped)) 
    Xy = csr_matrix(Xy)
   
    print("Final shape of Xy:", Xy.shape)

    return Xy
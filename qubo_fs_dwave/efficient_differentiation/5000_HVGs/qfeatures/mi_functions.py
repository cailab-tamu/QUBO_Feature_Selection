import numpy as np
from scipy.sparse import csr_matrix, issparse
from joblib import Parallel, delayed
import time

def binning_method(data, method="rice"):
    """
    Calculates the optimal binning for the given data based on the selected method.
    
    Parameters:
        data (array-like): The data for which to calculate the bins.
        method (str): The binning method. Options:
                      "auto", "square_root", "rice", "logarithmic", "freedman-diaconis", "scott".
    
    Returns:
        bins (ndarray or int): The bin edges (for "auto" method) or the number of bins (for other methods).
    """
    # Ensure data is a numpy array
    data = np.array(data)
    
    # Get the number of data points
    n = len(data)
    
    # If method is 'rice'
    if method == "rice":
        # Rice Rule: Number of bins = 2 * n^(1/3)
        return int(np.ceil(2 * n**(1/3)))

    # If method is 'logarithmic'
    elif method == "logarithmic":
        # Logarithmic Binning: Number of bins = log2(n) + 1
        return int(np.ceil(np.log2(n) + 1))
    
    # If method is 'square_root'
    elif method == "square_root":
        # Square Root Rule: Number of bins = sqrt(n)
        return int(np.ceil(np.sqrt(n)))
    
    # If method is 'freedman-diaconis'
    elif method == "freedman-diaconis":
        # Freedman-Diaconis Rule: Bin width = (2 * IQR) / n^(1/3)
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        bin_width = 2 * iqr / (n ** (1/3))
        return int(np.ceil((np.max(data) - np.min(data)) / bin_width))
    
    # If method is 'scott'
    elif method == "scott":
        # Scott's Rule: Bin width = (3.5 * std(data)) / (n^(1/3))
        bin_width_scott = (3.5 * np.std(data)) / (n ** (1/3))
        return int(np.ceil((np.max(data) - np.min(data)) / bin_width_scott))
    
    else:
        raise ValueError("Invalid method. Choose from: 'auto', 'square_root', 'rice', 'logarithmic', 'freedman-diaconis', 'scott'.")
    

def mutual_information_matrix(matrix, bins_method='rice', n_jobs=-1):
    """
    Computes the mutual information matrix in parallel, working directly with sparse matrices,
    and computes the full matrix (including the diagonal elements).
    """
    if not issparse(matrix):
        matrix = csr_matrix(matrix)

    n_features = matrix.shape[0]
    mi_matrix = np.zeros((n_features, n_features))

    def compute_pairwise_mi(i, j, matrix, bins_method='rice'):
        """
        Computes mutual information between row i and row j of the sparse matrix.
        """
        vi = matrix[i, :].toarray().flatten()
        vj = matrix[j, :].toarray().flatten()

        # Select binning method
        x_bins = binning_method(vi, bins_method)
        y_bins = binning_method(vj, bins_method)
        # Print the bins
        print(f"Bins for feature {i} (x):", x_bins)
        print(f"Bins for feature {j} (y):", y_bins)
                
        joint_counts, _, _ = np.histogram2d(vi, vj, bins=[x_bins, y_bins])
        ncounts = joint_counts.sum()
        if ncounts == 0:
            return 0  # No mutual information if no overlap
        joint_prob = joint_counts / ncounts + 1e-10

        marginal_i = joint_prob.sum(axis=1) + 1e-10
        marginal_j = joint_prob.sum(axis=0) + 1e-10
        
        joint_prob = joint_prob.flatten()
        h_xy =  -np.sum(joint_prob * np.log2(joint_prob))
        h_x =  -np.sum(marginal_i * np.log2(marginal_i))
        h_y =  -np.sum(marginal_j * np.log2(marginal_j))

        return float(h_x + h_y - h_xy)

    start_time = time.time()

    # Parallelizing the pairwise mutual information computation
    jobs = [(i, j) for i in range(n_features) for j in range(i+1, n_features)]  # Includes diagonal
    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_pairwise_mi)(i, j, matrix, bins_method='rice') for i, j in jobs
    )

    # Fill the matrix with the results
    for idx, (i, j) in enumerate(jobs):
        mi_matrix[i, j] = results[idx]
        mi_matrix[j, i] = results[idx]  # Exploit symmetry to avoid duplicate computation
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time for MI construction: {elapsed_time:.4f} seconds") 
    
    return mi_matrix
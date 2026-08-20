import numpy as np
from scipy.optimize import root_scalar, minimize
import time
from dwave.samplers import SteepestDescentSolver
from dwave.system import LeapHybridSampler  # Import for the Leap Hybrid Sampler
import dimod  # Import for working with QUBOs

# Convert QUBO matrix to dictionary format
def qubo_matrix_to_dict(Q):
    qubo = {}
    for i in range(Q.shape[0]):
        for j in range(Q.shape[1]):
            if Q[i, j] != 0:
                qubo[(i, j)] = Q[i, j]
    return qubo

# Define the optimize_qubo function (using SteepestDescentSolver to solve QUBO)
def optimize_qubo(alpha, R, J):
    """
    Computes the QUBO matrix and solves it using the Steepest Descent Solver.
    
    Parameters:
    - alpha (float): Mixing optimal parameter.
    - R (numpy.ndarray): Redundancy matrix (MI of feature-feature).
    - J (numpy.ndarray): Importance vector (MI of target-features).
    
    Returns:
    - n (int): Number of selected features (non-zero in solution vector).
    - result (dimod.SampleSet): QUBO solver results.
    """
    # Compute Q matrix
    Qmat = (1 - alpha) * R - alpha * np.diag(J)
    
    # Convert Q matrix to QUBO dictionary
    qubo = qubo_matrix_to_dict(Qmat)
    
    # Initialize the Steepest Descent Solver
    sampler = SteepestDescentSolver()
    
    # Solve the QUBO problem
    result = sampler.sample_qubo(qubo, num_reads=100)

    # Number of selected features (non-zero elements in the best solution vector)
    n = sum(result.first.sample.values())

    print(f"Solution's energy {result.first.energy} with alpha {alpha} and {n} features")  # Print energy of the solution

    return n, result


def annealing_process(MI_mat, K=100):

    # Redundancy matrix and importance vector from MI_mat and K
    R = MI_mat[:-1, :-1] / (K - 1)
    J = MI_mat[-1, :-1]

    # Start timing
    start_time = time.time()

    # Objective function to optimize
    def objective(alpha):
        return optimize_qubo(alpha, R, J)[0] - K

    # Try to find the optimal alpha using root_scalar
    try:
        root_method = 'brentq'  # Use 'brentq' | 'bisect' | 'secant' method for root_scalar
        result = root_scalar(objective, bracket=[0, 1], x0=0.5, method='bisect')
        alphasol = result.root
        print(f"Optimal alpha value from root_scalar: {alphasol}")

        _, xsol = optimize_qubo(alphasol, R, J)
        print(f"Energy from the solution: {xsol.first.energy}")  # Print energy of the solution
    except:
        # If root_scalar fails, fall back to minimize
        result = minimize(objective, x0=1, method='Nelder-Mead')
        alphasol = result.x[0]
        print(f"Optimal alpha value from minimize: {alphasol}")

        _, xsol = optimize_qubo(alphasol, R, J)
        print(f"Energy from the solution: {xsol.first.energy}")  # Print energy of the solution

    # End timing
    time_zerof = time.time() - start_time

    # Print computation time
    print(f"Time taken: {time_zerof:.2f} seconds")

    return alphasol, xsol


def quantum_test(alpha, MI_mat, df_features, K=100, mode='sa'):
    """
    Perform a quantum or simulated annealing test.

    Parameters:
        alpha (float): Trade-off parameter between redundancy and importance.
        MI_mat (numpy.ndarray): Mutual information matrix, with the last row/column being importance scores.
        df_features (pd.DataFrame): DataFrame containing feature information.
        K (int): Number of features (default=100).
        mode (str): Solver mode ('sa' for simulated annealing, 'qa' for quantum annealing).

    Returns:
        tuple: (new_df, qa_sol)
            new_df: Filtered and sorted DataFrame with selected features.
            qa_sol: Solution vector indicating selected features.
    """
    # Redundancy matrix and importance vector from MI_mat and K
    R = MI_mat[:-1, :-1] / (K - 1)
    J = MI_mat[-1, :-1]

    # Compute Q matrix
    Qmat = (1 - alpha) * R - alpha * np.diag(J)
    
    # Convert Q matrix to QUBO dictionary
    qubo = qubo_matrix_to_dict(Qmat)

    if mode == 'qa':    
        # Quantum (Hybrid) Solver
        print("Quantum (Hybrid) Solver Results...")

        # Using Leap Hybrid Sampler
        sampler_q = LeapHybridSampler()
        sampleset_q = sampler_q.sample_qubo(qubo)

        # Convert total runtime from microseconds (10^-6) to seconds
        runtime_all = sampleset_q.info.get('run_time', 0)
        runtime_all_seconds = runtime_all / 1_000_000

        # Convert QPU access time from microseconds (10^-6) to seconds
        annealing_time = sampleset_q.info.get('qpu_access_time', 0)
        annealing_time_seconds = annealing_time / 1_000_000

        print(f"Quantum Annealing time: {annealing_time_seconds} seconds")
        print(f"D-Wave Hybrid Solver time: {runtime_all_seconds} seconds")

    else:
        print("Simulated Annealing Solver Results...")
        # Using Steepest Descent Solver
        sampler_q = SteepestDescentSolver()
        sampleset_q = sampler_q.sample_qubo(qubo, num_reads=100)

    # Retrieve and analyze quantum results
    print("Energy: ", sampleset_q.first.energy)
    print("Occurrences: ", sampleset_q.first.num_occurrences)

    # Store selected features from quantum solution
    sample_dict_q = sampleset_q.first.sample

    # Reconstruct full-length decision vector matching original Qmat size
    n_features = Qmat.shape[0]
    qa_sol = np.zeros(n_features, dtype=int)
    for idx, val in sample_dict_q.items():
        qa_sol[idx] = int(val)

    # Multiply full selected vector by Qmat
    ener_per_feat = Qmat @ qa_sol

    # Add results to the DataFrame
    df_features = df_features.copy()
    df_features['feature_selected'] = qa_sol
    df_features['feature_score'] = ener_per_feat

    # Filter results
    filt_df_q = df_features[df_features['feature_selected'] > 0].copy()

    # Sort by 'feature_score'
    filt_df_q = filt_df_q.sort_values(by='feature_score', ascending=True).copy()

    new_df = filt_df_q.reset_index(drop=True)
    new_df.index = new_df.index.astype(str)

    return new_df, qa_sol

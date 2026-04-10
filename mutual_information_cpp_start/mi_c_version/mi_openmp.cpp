#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <omp.h>
#include "mmio.h"
#include <cstdlib>
#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/SparseCore> 


// Helper function for equal-width binning
std::vector<double> create_bin_edges(const std::vector<double>& vec, size_t nbins) {
    double min_value = *std::min_element(vec.begin(), vec.end());
    double max_value = *std::max_element(vec.begin(), vec.end());

    std::vector<double> bin_edges(nbins + 1);
    double bin_width = (max_value - min_value) / nbins;

    for (size_t i = 0; i <= nbins; ++i) {
        bin_edges[i] = min_value + i * bin_width;
    }
    bin_edges[nbins] = max_value;  // Ensure last bin includes max value
    return bin_edges;
}

// Compute the bin index for a given value (Python-like behavior)
inline int compute_bin_index(double value, const std::vector<double>& bin_edges) {
    size_t nbins = bin_edges.size() - 1;
    if (value == bin_edges[nbins]) {
        return static_cast<int>(nbins - 1);  // Assign max value to the last bin
    }
    auto it = std::upper_bound(bin_edges.begin(), bin_edges.end() - 1, value);
    return static_cast<int>(it - bin_edges.begin() - 1);
}

// Compute mutual information between two vectors
inline double compute_pairwise_mi(const std::vector<double>& vec1, const std::vector<double>& vec2, 
                                  const std::vector<double>& bin_edges1, const std::vector<double>& bin_edges2, 
                                  size_t nbins) {
    std::vector<std::vector<int>> joint_counts(nbins, std::vector<int>(nbins, 0));

    // Create joint histogram
    for (size_t idx = 0; idx < vec1.size(); ++idx) {
        int bin1 = compute_bin_index(vec1[idx], bin_edges1);
        int bin2 = compute_bin_index(vec2[idx], bin_edges2);
        joint_counts[bin1][bin2]++;
    }

    int total_count = vec1.size();

    // Normalize joint probabilities
    std::vector<std::vector<double>> joint_prob(nbins, std::vector<double>(nbins, 0.0));
    for (size_t i = 0; i < nbins; ++i) {
        for (size_t j = 0; j < nbins; ++j) {
            joint_prob[i][j] = static_cast<double>(joint_counts[i][j]) / total_count;
        }
    }

    // Compute marginal probabilities
    std::vector<double> marginal1(nbins, 0.0), marginal2(nbins, 0.0);
    for (size_t i = 0; i < nbins; ++i) {
        for (size_t j = 0; j < nbins; ++j) {
            marginal1[i] += joint_prob[i][j];
            marginal2[j] += joint_prob[i][j];
        }
    }

    // Compute mutual information (numerical stability corrections applied)
    double mi = 0.0;
    const double eps = 1e-12;  // Small epsilon to avoid log(0)

    for (size_t i = 0; i < nbins; ++i) {
        for (size_t j = 0; j < nbins; ++j) {
            if (joint_prob[i][j] > eps) {
                mi += joint_prob[i][j] * 
                      (std::log2(joint_prob[i][j]) - std::log2(marginal1[i] + eps) - std::log2(marginal2[j] + eps));
            }
        }
    }
    return mi;
}

// Compute error between Python lower triangular and C++ upper triangular matrices
double computeError(const Eigen::MatrixXd& cppMatrix, const Eigen::MatrixXd& pythonMatrix) {
    if (cppMatrix.rows() != pythonMatrix.rows() || cppMatrix.cols() != pythonMatrix.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for error computation.");
    }

    double error = 0.0;
    int count = 0;

    // Compare Python lower triangular (excluding diagonal) with C++ upper triangular
    for (int i = 1; i < cppMatrix.rows(); ++i) { // Start from row 1 to skip the diagonal
        for (int j = 0; j < i; ++j) { // Loop over the lower triangle in Python and upper triangle in C++
            error += std::pow(pythonMatrix(i, j) - cppMatrix(j, i), 2); // Compare Python lower (i, j) with C++ upper (j, i)
            count++;
        }
    }

    return std::sqrt(error / count); // RMSE
}

// Load sparse matrix from Matrix Market file
Eigen::SparseMatrix<double> loadMatrixMarketFile(const std::string& filename) {
    Eigen::SparseMatrix<double> matrix;
    FILE* file = fopen(filename.c_str(), "r");
    if (!file) {
        throw std::runtime_error("Failed to open Matrix Market file: " + filename);
    }

    MM_typecode matcode;
    if (mm_read_banner(file, &matcode) != 0) {
        throw std::runtime_error("Could not read Matrix Market banner.");
    }

    int rows, cols, nonzeros;
    if (mm_read_mtx_crd_size(file, &rows, &cols, &nonzeros) != 0) {
        throw std::runtime_error("Could not read matrix dimensions from Matrix Market file.");
    }

    std::vector<int> row_indices(nonzeros);
    std::vector<int> col_indices(nonzeros);
    std::vector<double> values(nonzeros);

    for (int i = 0; i < nonzeros; ++i) {
        // Check the return value of fscanf to ensure reading was successful
        if (fscanf(file, "%d %d %lf\n", &row_indices[i], &col_indices[i], &values[i]) != 3) {
            fclose(file);
            throw std::runtime_error("Error reading matrix element " + std::to_string(i));
        }

        row_indices[i]--;  // Adjust for 0-based indexing
        col_indices[i]--;
    }

    fclose(file);  // Close the file after reading

    matrix.resize(rows, cols);
    for (int i = 0; i < nonzeros; ++i) {
        matrix.insert(row_indices[i], col_indices[i]) = values[i];
    }

    return matrix;
}


// Main function with OpenMP parallelization
int main(int argc, char* argv[]) {
    Eigen::SparseMatrix<double> sparseMatrix;
    Eigen::SparseMatrix<double> pythonMI;

    int max_proc =  omp_get_num_procs();

    // Check if the number of threads is provided as a command line argument
    int num_threads = max_proc;  // Default number of threads
    if (argc > 1) {
        num_threads = std::atoi(argv[1]);  // Convert argument to integer
    }
    std::cout << "Requested procs : " << num_threads << std::endl;
    std::cout << "Max procs available: " << max_proc << std::endl;

    omp_set_num_threads(num_threads);  // Set the number of threads based on input

    try {
        sparseMatrix = loadMatrixMarketFile("sparse_matrix.mtx");
        pythonMI = loadMatrixMarketFile("sparse_matrix_mi.mtx");
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    size_t nbins = 20;
    Eigen::MatrixXd cppMI(sparseMatrix.rows(), sparseMatrix.rows());

    // Start measuring time for the parallel for loop
    auto start_time = std::chrono::high_resolution_clock::now();

    // Parallelize the outer loop using OpenMP
    #pragma omp parallel for schedule(dynamic, 10)
    for (int i = 0; i < sparseMatrix.rows(); ++i) {
        Eigen::VectorXd row_i = sparseMatrix.row(i).toDense();

        for (int j = i; j < sparseMatrix.rows(); ++j) {
            Eigen::VectorXd row_j = sparseMatrix.row(j).toDense();

            std::vector<double> vec1(row_i.data(), row_i.data() + row_i.size());
            std::vector<double> vec2(row_j.data(), row_j.data() + row_j.size());

            auto bin_edges1 = create_bin_edges(vec1, nbins);
            auto bin_edges2 = create_bin_edges(vec2, nbins);

            double mi_value = compute_pairwise_mi(vec1, vec2, bin_edges1, bin_edges2, nbins);

            cppMI(i, j) = mi_value;
            cppMI(j, i) = mi_value;
        }
    }

    // End measuring time for the parallel for loop
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;

    std::cout << "Time taken for the parallel for loop: " << elapsed_time.count() << " seconds" << std::endl;

    try {
        double error = computeError(cppMI, pythonMI);
        std::cout << "Error between C++ and Python MI matrices: " << error << std::endl;
    } catch (const std::invalid_argument& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }

    return 0;
}


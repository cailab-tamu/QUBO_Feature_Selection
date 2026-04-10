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
std::vector<double> create_bin_edges(const std::vector<double>& sorted_vec, size_t nbins) {
    std::vector<double> bin_edges;
    double min_value = sorted_vec.front();
    double max_value = sorted_vec.back();
    
    for (size_t i = 0; i <= nbins; ++i) {
        double bin_edge = min_value + (max_value - min_value) * i / nbins;
        bin_edges.push_back(bin_edge);
    }
    return bin_edges;
}

// Quantile-based binning (equal frequency)
std::vector<double> create_quantile_bin_edges(const std::vector<double>& sorted_vec, size_t nbins) {
    std::vector<double> bin_edges;
    size_t n = sorted_vec.size();
    
    for (size_t i = 0; i <= nbins; ++i) {
        size_t index = static_cast<size_t>(std::round(i * (n - 1) / static_cast<double>(nbins)));
        bin_edges.push_back(sorted_vec[index]);
    }
    return bin_edges;
}

// Logarithmic binning for wide range data
std::vector<double> create_log_bin_edges(const std::vector<double>& sorted_vec, size_t nbins) {
    std::vector<double> bin_edges;
    double min_value = sorted_vec.front();
    double max_value = sorted_vec.back();
    
    for (size_t i = 0; i <= nbins; ++i) {
        double log_min = std::log(min_value);
        double log_max = std::log(max_value);
        double bin_edge = std::exp(log_min + (log_max - log_min) * i / nbins);
        bin_edges.push_back(bin_edge);
    }
    return bin_edges;
}

// Dynamic binning strategy: Here, we dynamically decide the bin edges based on the data distribution
std::vector<double> create_dynamic_bin_edges(const std::vector<double>& sorted_vec, size_t nbins) {
    std::vector<double> bin_edges;
    size_t n = sorted_vec.size();

    // Calculate density-based binning (using interquartile range)
    double q1 = sorted_vec[n / 4];
    double q3 = sorted_vec[3 * n / 4];
    double iqr = q3 - q1; // Interquartile range (used for binning sensitivity)
    
    // The bins are distributed based on IQR and data spread
    double min_value = sorted_vec.front();
    double max_value = sorted_vec.back();

    for (size_t i = 0; i <= nbins; ++i) {
        double bin_edge = min_value + (max_value - min_value) * i / nbins;
        if (i > 0 && bin_edge < (q1 - 1.5 * iqr)) {
            bin_edge = (q1 - 1.5 * iqr); // Adjust for density-based bin edges
        }
        bin_edges.push_back(bin_edge);
    }
    return bin_edges;
}

// Compute the bin index for a given value
int compute_bin_index(double value, const std::vector<double>& bin_edges) {
    if (value < bin_edges.front()) {
        return -1;
    }
    if (value >= bin_edges.back()) {
        return bin_edges.size() - 2;
    }
    for (size_t i = 0; i < bin_edges.size() - 1; ++i) {
        if (value >= bin_edges[i] && value < bin_edges[i + 1]) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// Compute mutual information between two vectors
double compute_pairwise_mi(const std::vector<double>& vec1, const std::vector<double>& vec2, size_t nbins) {
    if (vec1.empty() || vec2.empty() || nbins == 0) {
        throw std::invalid_argument("Input vectors must not be empty, and nbins must be greater than 0.");
    }
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    std::vector<double> sorted_vec1 = vec1, sorted_vec2 = vec2;
    std::sort(sorted_vec1.begin(), sorted_vec1.end());
    std::sort(sorted_vec2.begin(), sorted_vec2.end());

    // Create bin edges based on the dynamic binning strategy
    std::vector<double> bin_edges1 = create_dynamic_bin_edges(sorted_vec1, nbins);
    std::vector<double> bin_edges2 = create_dynamic_bin_edges(sorted_vec2, nbins);

    // Create joint histogram
    std::vector<std::vector<int>> joint_counts(nbins, std::vector<int>(nbins, 0));

    for (size_t idx = 0; idx < vec1.size(); ++idx) {
        int bin1 = compute_bin_index(vec1[idx], bin_edges1);
        int bin2 = compute_bin_index(vec2[idx], bin_edges2);

        if (bin1 >= 0 && bin1 < nbins && bin2 >= 0 && bin2 < nbins) {
            joint_counts[bin1][bin2]++;
        }
    }

    // Calculate total count directly
    int total_count = std::accumulate(joint_counts.begin(), joint_counts.end(), 0,
        [](int sum, const std::vector<int>& row) {
            return sum + std::accumulate(row.begin(), row.end(), 0);
        });

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

    // Define epsilon for numerical stability
    const double eps0 = 1e-10;

    // Compute the joint entropy H(X, Y)
    double joint_entropy = 0.0;
    for (size_t i = 0; i < nbins; ++i) {
        for (size_t j = 0; j < nbins; ++j) {
            double joint_val = joint_prob[i][j] + eps0;
            joint_entropy -= joint_val * std::log2(joint_val);
        }
    }

    // Compute the entropy H(X)
    double entropy_x = 0.0;
    for (size_t i = 0; i < nbins; ++i) {
        double marginal_val_x = marginal1[i] + eps0;
        entropy_x -= marginal_val_x * std::log2(marginal_val_x);
    }

    // Compute the entropy H(Y)
    double entropy_y = 0.0;
    for (size_t j = 0; j < nbins; ++j) {
        double marginal_val_y = marginal2[j] + eps0;
        entropy_y -= marginal_val_y * std::log2(marginal_val_y);
    }

    // Compute the mutual information MI
    double mi = entropy_x + entropy_y - joint_entropy;

    return mi;
}

// Compute error between two matrices (lower triangle only)
double computeError(const Eigen::MatrixXd& matrix1, const Eigen::MatrixXd& matrix2) {
    if (matrix1.rows() != matrix2.rows() || matrix1.cols() != matrix2.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for error computation.");
    }

    double error = 0.0;
    int count = 0;  // Keep track of the number of comparisons

    // Loop over the lower triangle part (excluding diagonal) of the matrix
    for (int i = 1; i < matrix1.rows(); ++i) {  // Start from row 1 to skip the diagonal
        for (int j = 0; j < i; ++j) {  // Loop over columns to the left of the diagonal
            // Compute squared error between the values of the lower triangle
            error += std::pow(matrix1(i, j) - matrix2(i, j), 2);
            count++;
        }
    }

    // Calculate the root mean squared error (RMSE)
    return std::sqrt(error / count); // RMSE
}


// Load sparse matrix from Matrix Market file
Eigen::SparseMatrix<double> loadMatrixMarketFile(const std::string &filename) {
    Eigen::SparseMatrix<double> matrix;
    std::cout << "Loading matrix from file: " << filename << std::endl;

    // Open the file
    FILE* file = fopen(filename.c_str(), "r");
    if (!file) {
        throw std::runtime_error("Failed to open Matrix Market file: " + filename);
    }

    // Read matrix market header
    MM_typecode matcode;
    if (mm_read_banner(file, &matcode) != 0) {
        throw std::runtime_error("Could not read Matrix Market banner.");
    }

    // Check if the file contains a sparse matrix
    if (!mm_is_sparse(matcode)) {
        fclose(file);
        throw std::runtime_error("Matrix Market file is not sparse.");
    }

    int rows, cols, nonzeros;
    if (mm_read_mtx_crd_size(file, &rows, &cols, &nonzeros) != 0) {
        fclose(file);
        throw std::runtime_error("Could not read matrix dimensions from Matrix Market file.");
    }

    // Allocate space for matrix data (non-zero entries)
    std::vector<int> row_indices(nonzeros);
    std::vector<int> col_indices(nonzeros);
    std::vector<double> values(nonzeros);

    // Read matrix entries
    for (int i = 0; i < nonzeros; ++i) {
        if (fscanf(file, "%d %d %lf\n", &row_indices[i], &col_indices[i], &values[i]) != 3) {
            fclose(file);
            throw std::runtime_error("Error reading matrix entries from Matrix Market file.");
        }
        // Convert to 0-based indexing
        row_indices[i]--;
        col_indices[i]--;
    }

    fclose(file);

    // Create the Eigen SparseMatrix from the read values
    matrix.resize(rows, cols);
    for (int i = 0; i < nonzeros; ++i) {
        matrix.insert(row_indices[i], col_indices[i]) = values[i];
    }

    std::cout << "Matrix loaded with size: " << matrix.rows() << "x" << matrix.cols() << std::endl;

    // Debugging: Print a few matrix entries
    int max_entries_to_print = 10; // Limit the number of entries printed
    int count = 0;
    for (int k = 0; k < matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it) {
            if (count >= max_entries_to_print) {
                break; // Stop printing after max_entries_to_print
            }
            //std::cout << "Matrix entry: (" << it.row() << ", " << it.col() << ") = " << it.value() << std::endl;
            count++;
        }
        if (count >= max_entries_to_print) {
            break; // Stop iterating after max_entries_to_print
        }
    }

    return matrix;
}

// Main function
int main() {
    Eigen::SparseMatrix<double> sparseMatrix;
    Eigen::SparseMatrix<double> pythonMI;

    try {
        // Load sparse input matrix from Matrix Market file
        sparseMatrix = loadMatrixMarketFile("sparse_matrix.mtx");

        // Load Python-computed MI matrix
        pythonMI = loadMatrixMarketFile("sparse_matrix_mi.mtx");

        std::cout << "Sparse Matrix size: " << sparseMatrix.rows() << "x" << sparseMatrix.cols() << std::endl;
        std::cout << "Python MI Matrix size: " << pythonMI.rows() << "x" << pythonMI.cols() << std::endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }

    size_t nbins = 20;
    Eigen::MatrixXd cppMI(sparseMatrix.rows(), sparseMatrix.rows());

    // Start measuring time for the parallel for loop
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < sparseMatrix.rows(); ++i) {
        // Convert row `i` of the sparse matrix to a dense vector
        Eigen::VectorXd row_i = sparseMatrix.row(i).toDense();

        for (int j = i; j < sparseMatrix.rows(); ++j) {
            // Convert row `j` of the sparse matrix to a dense vector
            Eigen::VectorXd row_j = sparseMatrix.row(j).toDense();

            // Transform Eigen::VectorXd to std::vector<double> for the MI computation
            std::vector<double> vec1(row_i.data(), row_i.data() + row_i.size());
            std::vector<double> vec2(row_j.data(), row_j.data() + row_j.size());

            // Compute the mutual information
            double mi_value = compute_pairwise_mi(vec1, vec2, nbins);

            // Store the MI value (symmetrical matrix)
            cppMI(i, j) = mi_value;
            cppMI(j, i) = mi_value;

            // Print the MI value for the pair (i, j)
            //std::cout << "MI between row " << i << " and row " << j << " = " << mi_value << std::endl;
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

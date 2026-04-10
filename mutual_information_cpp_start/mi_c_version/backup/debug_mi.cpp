#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Sparse>
#include <numeric>
#include <stdexcept>

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

int main() {
    // Define a 4x5 matrix with a few zeros
    Eigen::SparseMatrix<double> sparse_matrix(4, 5);
    sparse_matrix.insert(0, 0) = 1;
    sparse_matrix.insert(0, 1) = 2;
    sparse_matrix.insert(0, 2) = 3;
    sparse_matrix.insert(1, 0) = 4;
    sparse_matrix.insert(1, 3) = 5;
    sparse_matrix.insert(2, 2) = 6;
    sparse_matrix.insert(3, 4) = 7;
    sparse_matrix.makeCompressed();

    std::cout << "Sparse matrix:\n";
    for (int k = 0; k < sparse_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(sparse_matrix, k); it; ++it) {
            std::cout << "(" << it.row() << ", " << it.col() << ") = " << it.value() << "\n";
        }
    }

    // Extract two rows as vectors
    std::vector<double> vec1 = {1, 2, 3, 0, 0};  // Row 0
    std::vector<double> vec2 = {4, 0, 6, 0, 0};  // Row 1

    // Compute MI between the two vectors
    size_t nbins = 5; // Number of bins for quantiles
    try {
        double mi = compute_pairwise_mi(vec1, vec2, nbins);
        std::cout << "Mutual Information (MI) between row 0 and row 1: " << mi << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }

    return 0;
}

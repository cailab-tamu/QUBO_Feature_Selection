

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <mpi.h>

// Entropy calculation
double entropy(const std::vector<double>& prob) {
    double result = 0.0;
    for (double p : prob) {
        if (p > 0) {
            result -= p * std::log2(p);
        }
    }
    return result;
}

// Matrix chunk handling and mutual information computation
double compute_pairwise_mi(const std::vector<std::vector<double>>& matrix, int i, int j, int nbins, int rank, int size) {
    int num_rows = matrix.size();
    int rows_per_process = num_rows / size;
    int start_row = rank * rows_per_process;
    int end_row = (rank == size - 1) ? num_rows : start_row + rows_per_process;

    double min_vi = std::numeric_limits<double>::infinity();
    double max_vi = -std::numeric_limits<double>::infinity();
    for (int row = start_row; row < end_row; ++row) {
        min_vi = std::min(min_vi, matrix[row][i]);
        max_vi = std::max(max_vi, matrix[row][i]);
    }

    double min_vj = std::numeric_limits<double>::infinity();
    double max_vj = -std::numeric_limits<double>::infinity();
    for (int row = start_row; row < end_row; ++row) {
        min_vj = std::min(min_vj, matrix[row][j]);
        max_vj = std::max(max_vj, matrix[row][j]);
    }

    std::vector<std::vector<int>> joint_counts(nbins, std::vector<int>(nbins, 0));

    for (int row = start_row; row < end_row; ++row) {
        double scaled_i = (matrix[row][i] - min_vi) / (max_vi - min_vi);
        double scaled_j = (matrix[row][j] - min_vj) / (max_vj - min_vj);

        int bin_i = std::min(nbins - 1, static_cast<int>(scaled_i * nbins));
        int bin_j = std::min(nbins - 1, static_cast<int>(scaled_j * nbins));
        joint_counts[bin_i][bin_j]++;
    }

    std::vector<std::vector<int>> global_joint_counts(nbins, std::vector<int>(nbins, 0));
    MPI_Reduce(&joint_counts[0][0], &global_joint_counts[0][0], nbins * nbins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double joint_sum = 0;
        for (const auto& row : global_joint_counts) {
            joint_sum += std::accumulate(row.begin(), row.end(), 0);
        }

        if (joint_sum == 0) {
            return 0.0;
        }

        std::vector<std::vector<double>> joint_prob(nbins, std::vector<double>(nbins, 0.0));
        for (int m = 0; m < nbins; ++m) {
            for (int n = 0; n < nbins; ++n) {
                joint_prob[m][n] = global_joint_counts[m][n] / joint_sum;
            }
        }

        std::vector<double> marginal_i(nbins, 0.0);
        std::vector<double> marginal_j(nbins, 0.0);
        for (int m = 0; m < nbins; ++m) {
            for (int n = 0; n < nbins; ++n) {
                marginal_i[m] += joint_prob[m][n];
                marginal_j[n] += joint_prob[m][n];
            }
        }

        double h_xy = entropy(joint_prob[0]);
        double h_x = entropy(marginal_i);
        double h_y = entropy(marginal_j);

        return h_x + h_y - h_xy;
    }

    return 0.0;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::vector<double>> matrix = {
        {0.1, 0.2, 0.3, 0.4},
        {0.5, 0.6, 0.7, 0.8},
        {0.2, 0.3, 0.4, 0.5},
        {0.3, 0.4, 0.5, 0.6},
        {0.7, 0.8, 0.9, 1.0}
    };

    int i = 0, j = 1;
    int nbins = 5;

    double mi = compute_pairwise_mi(matrix, i, j, nbins, rank, size);

    double global_mi;
    MPI_Reduce(&mi, &global_mi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Mutual Information: " << global_mi << std::endl;
    }

    MPI_Finalize();
    return 0;
}

% Load necessary data
load Data_hESC_EC_day1_5000g.mat
load Tqubo__R0_f50monocle3_pseudotime_HVG_remonocle_5000_5k_cells.mat
load R0.mat

% Set parameters
K = 20;
genes = sce.g;
alphasol = Tqubo.alphasol;
R = R0(1:end-1, 1:end-1) / (K - 1);
J = R0(end, 1:end-1);

% Recompute the Q matrix (balance of R and J) using trained alpha
Q = (1 - alphasol) * R - alphasol * diag(J);

% Generate the idxg based on Tqubo.selectedGenes and genes
idxg = zeros(K, 1);
for ig = 1:length(idxg)
    idxg(ig) = find(genes == Tqubo.selectedGenes(ig)); % Map to corresponding index
end

% Extract Qsub matrix (before shuffling)
Qsub_sorted = Q(idxg, idxg);  % Sorted Qsub matrix (based on original idxg)

% Get all possible permutations of x with 3 ones in a vector of length 20
combinations = nchoosek(1:20, 3);

% Initialize an array to store the energy values
E_icomb = zeros(size(combinations, 1), 1);

% Loop through all combinations to compute energy values
for i = 1:size(combinations, 1)
    % Create a vector x_icomb with 1s at the selected positions
    x_icomb = zeros(20, 1);
    x_icomb(combinations(i, :)) = 1;
    
    % Calculate the energy E_icomb = x_icomb' * (Qsub_sorted * x_icomb)
    Evec = Qsub_sorted * x_icomb;
    E_icomb(i) = x_icomb' * Evec;
end

% Reverse the energy values (for sorted case)
E_icomb_reversed = flip(E_icomb);

% Plot the reversed 1D energy values
figure;
hold on;
plot(E_icomb_reversed, '-o', 'MarkerSize', 4, 'DisplayName', 'Raw Energy Values (Reversed)');
smoothed_E_reversed = smooth(E_icomb_reversed, 0.1, 'loess');
plot(smoothed_E_reversed, '-', 'LineWidth', 2, 'DisplayName', 'Smoothed Curve (Reversed)');
xlabel('Combination Index');
ylabel('Energy Value');
title('Energy Landscape for All Combinations of 20 Genes and Selecting 3');
legend('show');
grid on;
hold off;

%% Shuffle Qsub matrix for 2D energy plot
% Shuffle idxg (random permutation of idxg)
idxg_shuffled = idxg(randperm(length(idxg))); % Shuffle idxg randomly

% Permute Q using shuffled idxg for both rows and columns
Qsub_shuffled = Q(idxg_shuffled, idxg_shuffled);  % Apply the random permutation to both rows and columns

% Initialize an array to store the energy values for shuffled case
E_icomb_shuffled = zeros(size(combinations, 1), 1);

% Loop through all combinations to compute energy values for shuffled Q
for i = 1:size(combinations, 1)
    % Create a vector x_icomb with 1s at the selected positions
    x_icomb = zeros(20, 1);
    x_icomb(combinations(i, :)) = 1;
    
    % Calculate the energy E_icomb = x_icomb' * (Qsub_shuffled * x_icomb)
    Evec = Qsub_shuffled * x_icomb;
    E_icomb_shuffled(i) = x_icomb' * Evec;
end

% Reverse the energy values (for shuffled case)
E_icomb_shuffled_reversed = flip(E_icomb_shuffled);

% Reshape the reversed energy values for 2D grid (ensure reshaped array fits grid dimensions)
n1 = 30; n2 = 38;  % Set grid dimensions
pZ = reshape(E_icomb_shuffled_reversed, n2, n1);  % Reshaping to fit (n2, n1) grid

% Create a 3D surface plot for the shuffled energy values
hx = gui.myFigure;
[pX, pY] = meshgrid(1:n1, 1:n2);  % Define grid for X and Y axes
s = surf(pX, pY, pZ, 'EdgeColor', 'none');  % Create the surface plot

% Set plot labels and view options
xlabel('Combination Index 1');
ylabel('Combination Index 2');
zlabel('Energy Value');
title('Shuffled Energy Plot of Combinations');
box on;
view(3);  % 3D view for better visualization

hx.show;  % Display the figure

% Load necessary data
load ../Data_hESC_EC_day1_5000g.mat
load ../selected_50f_results/Tqubo_f50_monocle3_pseudotime_HVG_5000_cells_4697.mat
load ../R0.mat

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

%qubogenes = Tqubo.selectedGenes(randperm(numel(Tqubo.selectedGenes)));
qubogenes = Tqubo.selectedGenes;
for ig = 1:K
    idxg(ig) = find(genes == qubogenes(ig)); % Map to corresponding index
end
rng("default")
idxg = idxg(randperm(numel(idxg)));

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
f = figure;
hold on;
plot(E_icomb_reversed, '-o', 'MarkerSize', 4, 'DisplayName', 'Raw Energy Values');
smoothed_E_reversed = smooth(E_icomb_reversed, 0.1, 'loess');
plot(smoothed_E_reversed, '-', 'LineWidth', 2, 'DisplayName', 'Smoothed Curve');
xlabel('Combination index','FontSize',14);
ylabel('Energy value','FontSize',14);
title({'Energy Landscape 1D', ' All combinations of 20 genes selecting 3'}, 'FontSize', 15);
lgd = legend('show', 'Location','best');
fontsize(lgd, 11,'points')
grid on;
hold off;

% Save as high-quality PNG using print
filename = 'raw_energy_combinations.svg';
print(f, filename, '-dsvg');
%% Plot sorted selection in 1D energy values

% Generate the idxg based on Tqubo.selectedGenes and genes
idxg = zeros(K, 1);

%qubogenes = Tqubo.selectedGenes(randperm(numel(Tqubo.selectedGenes)));
qubogenes = Tqubo.selectedGenes;
for ig = 1:K
    idxg(ig) = find(genes == qubogenes(ig)); % Map to corresponding index
end

% Extract Qsub matrix (before shuffling)
Qsub_sorted = Q(idxg, idxg);  % Sorted Qsub matrix (based on original idxg)

% Get all possible permutations of x with 3 ones in a vector of length 20
combinations = nchoosek(1:20, 3);

% Initialize an array to store the energy values
E_icomb2 = zeros(size(combinations, 1), 1);

% Loop through all combinations to compute energy values
for i = 1:size(combinations, 1)
    % Create a vector x_icomb with 1s at the selected positions
    x_icomb = zeros(20, 1);
    x_icomb(combinations(i, :)) = 1;
    
    % Calculate the energy E_icomb = x_icomb' * (Qsub_sorted * x_icomb)
    Evec = Qsub_sorted * x_icomb;
    E_icomb2(i) = x_icomb' * Evec;
end

% Reverse the energy values (for sorted case)
E_icomb_reversed2 = flip(E_icomb2);

% Plot the reversed 1D energy values
f = figure;
hold on;
plot(E_icomb_reversed2, '-o', 'MarkerSize', 4, 'DisplayName', 'Raw Energy Values');
smoothed_E_reversed = smooth(E_icomb_reversed2, 0.1, 'loess');
plot(smoothed_E_reversed, '-', 'LineWidth', 2, 'DisplayName', 'Smoothed Curve');
xlabel('Combination index','FontSize', 14);
ylabel('Energy value','FontSize', 14);
title({'Energy minimization process', ' All combinations of 20 genes selecting 3'}, 'FontSize', 15);
lgd = legend('show', 'Location','best');
fontsize(lgd, 11,'points')
grid on;
hold off;

% Save as high-quality PNG using print
filename = 'energy_combinations_min.svg';
print(f, filename, '-dsvg');

%% Shuffle Qsub matrix for 2D energy plot
%{ 
%Shuffle idxg (random permutation of idxg)
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

%}
% Reshape the reversed energy values for 2D grid (ensure reshaped array fits grid dimensions)
setid = 1; 
smoothed = false;
useinterp = true;

switch setid
    case 1
        n1 = 30; n2 = 38;  % Set grid dimensions
    case 2
        n1 = 20; n2 = 57;  % Set grid dimensions
    case 3
        n1 = 10; n2 = 114;  % Set grid dimensions
end

if smoothed
    pZ = reshape(smoothed_E_reversed, n2, n1);  % Reshaping to fit (n2, n1) grid
else
    pZ = reshape(E_icomb_reversed, n2, n1);  % Reshaping to fit (n2, n1) grid
end

% Create a 3D surface plot for the shuffled energy values
hx = gui.myFigure;
[pX, pY] = meshgrid(1:n1, 1:n2);  % Define grid for X and Y axes

if useinterp
    % Create a finer grid for interpolation
    [Xq, Yq] = meshgrid(1:0.1:n1, 1:0.1:n2); % Finer grid
    Zq = interp2(pX, pY, pZ, Xq, Yq, 'spline'); % Use 'spline' for smooth interpolation
    s = surf(Xq, Yq, Zq, 'EdgeColor', 'none','FaceAlpha', 0.5);
else
    s = surf(pX, pY, pZ, 'EdgeColor', 'none');  % Create the surface plot
end

%box on;
colormap jet
view(3);  % 3D view for better visualization
hx.show;  % Display the figure

ax = hx.FigHandle.Children; % Access the Axes object ( not working again)

xlabel('Combination index 1', 'FontSize', 14); % Use 'Parent' to specify axes
ylabel('Combination index 2', 'FontSize', 14); % Use 'Parent' to specify axes
zlabel('Energy value', 'FontSize', 14);
title({'Energy landscape 3D', ' All combinations of 20 genes selecting 3'}, 'FontSize', 15);

% Rotate the labels
xlabel_handle = ax.XLabel; % Access the XLabel object directly
ylabel_handle = ax.YLabel; % Access the YLabel object directly
set(xlabel_handle, 'Rotation', 16);
set(ylabel_handle, 'Rotation', -24);

% Save as high-quality PNG using print
filename_svg = 'energy_combinations_landscape.svg';
% Save the plot as an SVG file
print(ax, filename_svg, '-dsvg');
%% Energy path plot
% Energy landscape
load("../R0.mat")
load("../selected_500f_results/Tqubo_f500_monocle3_pseudotime_HVG_5000_cells_4697.mat")
load("../selected_500f_results/Tml_f500_monocle3_pseudotime_HVG_5000_cells_4697.mat")

K = length( Tqubo.selectedGenes); 
energy_path3(R0, Tqubo.selectedGenes, Tml{1}.selectedGenes, Tml{4}.selectedGenes, genes, K,...
                 Tqubo.alphasol, 'energy_path2_500.svg');

%% Not used...
% a=dec2bin([18 29 20;1 2 3]', 5)'-'0'
a = dec2bin(combinations', 5)'-'0';
a = a(:);
b = reshape(a, [15 length(a)/15])';
rng(500);
[eb] = tsne(b);
x = eb(:,1); y = eb(:,2);
z = E_icomb_reversed;
[xq, yq] = meshgrid(linspace(min(x), max(x), 50), linspace(min(y), max(y), 50)); 
zq = griddata(x, y, z, xq, yq, 'natural');



%{
    zq = fillmissing(zq, 'linear', 2); % Fill NaNs linearly along rows
    zq = fillmissing(zq, 'linear', 1); % Then fill along columns
    % zq(isnan(zq)) = min(zq(:)); % mean(zq(:), 'omitnan');
    [Xq, Yq] = meshgrid(1:0.1:size(xq, 1), 1:0.1:size(xq, 2)); % Finer grid
    [Xq, Yq] = meshgrid(1:0.1:size(xq, 1), 1:0.1:size(xq, 2)); % Finer grid
    Zq = gridfit(xq, yq, zq, 1:0.1:size(xq, 1), 1:0.1:size(xq, 2));
    % Zq = interp2(xq, yq, zq, Xq, Yq, 'spline'); % Use 'spline' for smooth interpolation    
    surf(Xq, Yq, Zq, 'EdgeColor', 'none');
%}

hx = gui.myFigure;
% surf(Xq, Yq, Zq); % Surface plot
surf(xq, yq, zq, 'EdgeColor', 'none');
shading interp;   % Smooth shading
xlabel('Combination Index 1');
ylabel('Combination Index 2');
zlabel('Energy Value');
title('Energy landscape for combinations of 3 out of 20 genes');
colormap(flipud(hot));
colorbar
view(3);  % 3D view for better visualization
hx.show;  % Display the figure


%{
%https://arxiv.org/abs/2411.14708
% colormap(flipud(jet));
% Rugged surface of a 5D Sphere function
when inputs are represented as Gemini embeddings of
dimension 6K+, post-processed by t-SNE into 2D space.

hx = gui.myFigure;
s = surf(z(:,1), z(:,2), E_icomb_reversed, 'EdgeColor', 'none');  % Create the surface plot
xlabel('Combination Index 1');
ylabel('Combination Index 2');
zlabel('Energy Value');
title('Energy landscape for combinations of 3 out of 20 genes');
colormap jet
view(3);  % 3D view for better visualization
hx.show;  % Display the figure
%}

% Define the Q matrix (20x20)
% Redundancy matrix (scaled)
load Data_hESC_EC_day1_5000g.mat
load Tqubo__R0_f50monocle3_pseudotime_HVG_remonocle_5000_5k_cells.mat
load R0.mat
K = 20; 
genes = sce.g;
alphasol = Tqubo.alphasol;
R = R0(1:end-1, 1:end-1) / (K - 1);

% Importance vector
J = R0(end, 1:end-1);

% Recompute Q matrix (balance of R and J) using trained alpha
Q = (1 - alphasol) * R - alphasol * diag(J);

idxg = zeros(K,1);
for ig = 1:length(idxg)
    idxg(ig) = find( genes == Tqubo.selectedGenes(ig) );
end

Qsub = Q(idxg,idxg);

% Get all possible permutations of x
% Since x only has 3 ones, we can use the combination function to find all positions for ones
combinations = nchoosek(1:20, 3);

% Initialize an array to store the energy values
E_icomb = zeros(size(combinations, 1), 1);

% Loop through all combinations
for i = 1:size(combinations, 1)
    % Create a new vector x_icomb with 1s at the selected positions
    x_icomb = zeros(20, 1);
    x_icomb(combinations(i, :)) = 1;
    
    % Calculate the energy E_icomb = sum(Q * x_icomb)
    Evec = Qsub*x_icomb;
    E_icomb(i) = x_icomb'*Evec;
end

% Reverse the energy values
E_icomb_reversed = flip(E_icomb);

% Plot the reversed energy values
figure;
hold on;  % Allow multiple plots on the same figure
plot(E_icomb_reversed, '-o', 'MarkerSize', 4, 'DisplayName', 'Raw Energy Values (Reversed)');  % Raw energy values (reversed)

% Smooth the reversed energy values
smoothed_E_reversed = smooth(E_icomb_reversed, 0.1, 'loess');  % Smooth the reversed values
plot(smoothed_E_reversed, '-', 'LineWidth', 2, 'DisplayName', 'Smoothed Curve (Reversed)');  % Smoothed curve (reversed)

% Add labels, legend, and grid
xlabel('Combination Index');
ylabel('Energy Value');
title('Energy Landscape for All Combinations of 20 Genes and Selecting 3');
legend('show');
grid on;
hold off;

%%
n1=30; n2=38;
n1=57; n2=20;
% n1=10; n2=114;

assert(~isempty(which('gui.myFigurex')),'needs scGEAToolbox')

hx = gui.myFigure;
[pX, pY] = meshgrid(1:n1, 1:n2);
% Z = reshape(randn(nchoosek(20, 3), 1), 38, 30);
pZ = reshape(smoothed_E_reversed, n2, n1);
s = surf(pX,pY,pZ,'EdgeColor','none');  %,'FaceAlpha',0.5);
% s(1).EdgeColor = 'none';
% sc(2).ZLocation = 'zmax';
xlabel('Combination Index 1');
ylabel('Combination Index 2');
zlabel('Energy Value');
box on
view(3)
hx.show;


%%
hx=gui.myFigure;
[pX, pY] = meshgrid(1:n1, 1:n2);
pZ=Z;
load earth
s=surface(pX,pY,pZ,X,'FaceColor','texturemap', ...
    'EdgeColor','none','CDataMapping','direct');
colormap(map)
view(3)
hx.show;

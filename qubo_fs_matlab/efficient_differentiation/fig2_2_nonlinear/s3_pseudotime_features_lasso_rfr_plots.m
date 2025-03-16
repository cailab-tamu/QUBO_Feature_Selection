load('../Data_hESC_EC_day1_5000g.mat')
g = sce.g;
X = full(sc_transform(sce.X, "type", "PearsonResiduals"));
K = 50;

%Tqubo = qfeatures_qubo_base(X, g, y, K, false);
g1 = Tqubo.selectedGenes;
%Tml = mlfeatures_base(X, g, y, K, 1);  
g2 = Tml{1}.selectedGenes;
% Random forest regression
g3 = Tml{4}.selectedGenes;

rfr_u_lasso = union(g2, g3);
qubo_only = setdiff(g1, rfr_u_lasso,'stable'); % This keeps g1 order
writematrix(qubo_only','qubo_only_features.txt')
qubo_only = g1; % Using all QUBO features 

% QUBO only features
idx_qubo = zeros(length(qubo_only), 1);
for i = 1:length(qubo_only)
    idx_qubo(i) = find(sce.g == qubo_only(i));
end
g_qubo = sce.g(idx_qubo);
X_qubo = X(idx_qubo,:);

% ML features
idx_ml = zeros(length(rfr_u_lasso), 1);
for i = 1:length(rfr_u_lasso)
    idx_ml(i) = find(sce.g == rfr_u_lasso(i));
end
g_ml = sce.g(idx_ml);
X_ml = X(idx_ml,:);

% Predictor
cell_type_target = 'splinefit_pseudotime';

% Preparing target predictor y from pseudo-time values per cell
idx = find(contains(sce.list_cell_attributes(1:2:end), cell_type_target));
if isempty(idx), return; end
t = sce.list_cell_attributes{idx*2};

% Smoothing parameter
sp = 0.75;

tic;
% QUBO fitting
ngene = size(X_qubo, 1);
ncell = size(X_qubo, 2);
y_qubo_fit = zeros(ncell, ngene);
t_qubo_sort = zeros(ncell, ngene);
for ig = 1:ngene
    [y_qubo_fit(1:ncell, ig), idx] = ...
                        loess_smoothing(t, X_qubo(ig,:)', sp);
    t_qubo_sort(1:ncell, ig) = t(idx);
end
toc;

tic;
% ML fitting
ngene = size(X_ml, 1);
ncell = size(X_ml, 2);
y_ml_fit = zeros(ncell, ngene);
t_ml_sort = zeros(ncell, ngene);
for ig = 1:ngene
    [y_ml_fit(1:ncell, ig), idx] = ...
                        loess_smoothing(t, X_ml(ig,:)', sp);
    t_ml_sort(1:ncell, ig) = t(idx);
end
toc;


%% Plot ML features in one figure
f_ml = figure;
hold on
num_ml_features = size(y_ml_fit, 2); % Get the number of ML features
for k = 1:num_ml_features
    plot(t_ml_sort(:,k), y_ml_fit(:,k),'LineWidth',2,'Color', 'k'); % Using blue for ML in its own figure
end
xlim([0 max(t)]);
ylim([-2 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('ML Features')
hold off;

%% Plot both ML and QUBO features in a different figure
f_combined = figure;
hold on
num_ml_features_combined = size(y_ml_fit, 2);
for k = 1:num_ml_features_combined
    plot(t_ml_sort(:,k), y_ml_fit(:,k),'LineWidth',2,'Color', 'k');
end
num_qubo_features = size(y_qubo_fit, 2);
for k = 1:num_qubo_features
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',0.1,'Color', 'g');
end
xlim([0 max(t)]);
ylim([-2 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('ML (Black) and QUBO (Green) Features')
hold off;


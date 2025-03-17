load('../GSE134839_cleandata_lite.mat')
g = sce.g;
X = sc_transform(sce.X);
K = 50;

load('../selected_100f_results/Tqubo_f100_manual_pseudotime_Genes_9661_cells_1422_cells_1422.mat')
g_qubo_sel = Tqubo.selectedGenes(1:K);

load('../selected_100f_results/Tml_f100_manual_pseudotime_Genes_9661_cells_1422.mat')
% LASSO features
g_lasso_sel = Tml{1}.selectedGenes(1:K);
% Random forest regression features
g_rfr_sel = Tml{4}.selectedGenes(1:K);

rfr_u_lasso = union(g_lasso_sel, g_rfr_sel);
qubo_only = setdiff(g_qubo_sel, rfr_u_lasso,'stable'); % This keeps g_qubo_sel order

% QUBO only features
idx_qubo_only = zeros(length(qubo_only), 1);
for i = 1:length(qubo_only)
    idx_qubo_only(i) = find(sce.g == qubo_only(i));
end

% QUBO features
idx_qubo = zeros(length(g_qubo_sel), 1);
for i = 1:length(g_qubo_sel)
    idx_qubo(i) = find(sce.g == g_qubo_sel(i));
end
g_qubo = sce.g(idx_qubo);
X_qubo = X(idx_qubo,:);

% LASSO features
idx_lasso = zeros(length(g_lasso_sel), 1);
for i = 1:length(g_lasso_sel)
    idx_lasso(i) = find(sce.g == g_lasso_sel(i));
end
g_lasso = sce.g(idx_lasso);
X_lasso = X(idx_lasso,:);

% RFR features
idx_rfr = zeros(length(g_rfr_sel), 1);
for i = 1:length(g_rfr_sel)
    idx_rfr(i) = find(sce.g == g_rfr_sel(i));
end
g_rfr = sce.g(idx_rfr);
X_rfr = X(idx_rfr,:);

% Predictor
cell_type_target = 'manual_pseudotime';

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
% LASSO fitting
ngene = size(X_lasso, 1);
ncell = size(X_lasso, 2);
y_lasso_fit = zeros(ncell, ngene);
t_lasso_sort = zeros(ncell, ngene);
for ig = 1:ngene
    [y_lasso_fit(1:ncell, ig), idx] = ...
                        loess_smoothing(t, X_lasso(ig,:)', sp);
    t_lasso_sort(1:ncell, ig) = t(idx);
end
toc;

tic;
% RFR fitting
ngene = size(X_rfr, 1);
ncell = size(X_rfr, 2);
y_rfr_fit = zeros(ncell, ngene);
t_rfr_sort = zeros(ncell, ngene);
for ig = 1:ngene
    [y_rfr_fit(1:ncell, ig), idx] = ...
                        loess_smoothing(t, X_rfr(ig,:)', sp);
    t_rfr_sort(1:ncell, ig) = t(idx);
end
toc;

% Descide min and max
rfr_max = max(max(y_rfr_fit));
lasso_max = max(max(y_lasso_fit));
qubo_max = max(max(y_qubo_fit));
ymax = round(max(rfr_max, max(lasso_max, qubo_max)));

rfr_min = min(min(y_rfr_fit));
lasso_min = min(min(y_lasso_fit));
qubo_min = min(min(y_qubo_fit));
ymin = round(min(rfr_min, min(lasso_min, qubo_min)));

%% Plot LASSO features in one figure
f = figure;
hold on
ng = size(y_lasso_fit, 2); % Get the number of ML features
for k = 1:ng
    plot(t_lasso_sort(:,k), y_lasso_fit(:,k),'LineWidth',2,'Color', 'k'); % Using blue for ML in its own figure
end
xlim([0 max(t)]);
ylim([ymin ymax]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('LASSO Features')
hold off;
f.Position(3) = 545;
f.Position(4) = 266;
filename = 'lasso_features.png';

% Save as high-quality PNG using print
resolution = 300; 
print(f, filename, '-dpng', sprintf('-r%d', resolution));

%% Plot RFR features in one figure
f = figure;
hold on
ng = size(y_rfr_fit, 2); % Get the number of ML features
for k = 1:ng
    plot(t_rfr_sort(:,k), y_rfr_fit(:,k),'LineWidth',2,'Color', 'k'); % Using blue for ML in its own figure
end
xlim([0 max(t)]);
ylim([ymin ymax]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('RFR Features')
hold off;
f.Position(3) = 545;
f.Position(4) = 266;
filename = 'rfr_features.png';

% Save as high-quality PNG using print
resolution = 300; 
print(f, filename, '-dpng', sprintf('-r%d', resolution));

%% Plot QUBO unique and all features
f = figure;
hold on
ng = size(y_qubo_fit, 2); % Get the number of ML features
for k = 1:ng
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',2,'Color', 'k'); % Using blue for ML in its own figure
end

idx_qubo_only = ismember(g_qubo_sel, qubo_only);
t_qubo_unique_sort = t_qubo_sort(:,idx_qubo_only);
y_qubo_unique_fit = y_qubo_fit(:,idx_qubo_only);

ng = size(y_qubo_unique_fit, 2);
for k = 1:ng
    plot(t_qubo_unique_sort(:,k), y_qubo_unique_fit(:,k),'LineWidth',0.5,'Color', 'g');
end
xlim([0 max(t)]);
ylim([ymin ymax]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('QUBO Features')
hold off;
f.Position(3) = 545;
f.Position(4) = 266;
filename = 'qubo_features.png';

% Save as high-quality PNG using print
resolution = 300; 
print(f, filename, '-dpng', sprintf('-r%d', resolution));

%% Unique QUBO genes with feature names 
my_genes = qubo_only;
% Create a 6x5 figure
f = figure;
tiledlayout(3, 5); % Use tiledlayout for better subplot management

for idx = 1:length(my_genes)
    % subplot(2, 5, i); % Create subplot
    nexttile;
    hold on;
    % Plot all QUBO genes in gray
    for k = 1:50
        plot(t_qubo_sort(:,k), y_qubo_fit(:,k), ...
            'LineWidth', 2, 'Color', 'k'); %[0.35, 0.35, 0.35]);
        %plot(t_qubo_sort(:,k), y_qubo_fit(:,k), 'LineWidth', 2, 'Color', 'k');
    end
    % Plot the selected gene in red
    gene_index = find(strcmp(g_qubo, my_genes(idx)));
    plot(t_qubo_sort(:,gene_index), ...
        y_qubo_fit(:,gene_index), 'LineWidth', 2, ...
        'Color', 'r');
    xlim([-0.06 max(t)+0.04]);
    ylim([ymin ymax]);
    
    %xlabel('Pseudotime')
    %ylabel('Standardized Expression')
    box on
    title(sprintf('%s', my_genes(idx)));
end
f.Position(3:4) = [1364 480];

% Save as SVG
filename = 'qubo_unique_genes_figure.png';

% Save as high-quality PNG using print
resolution = 300; 
print(f, filename, '-dpng', sprintf('-r%d', resolution));
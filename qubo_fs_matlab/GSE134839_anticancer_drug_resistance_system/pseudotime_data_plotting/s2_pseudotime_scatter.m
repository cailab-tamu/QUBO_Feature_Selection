tic;
load('../GSE134839_cleandata_lite.mat')

infile = 'genes_all_lasso_manual_pseudo.txt'; 
highlightg = false; ttxt = 'LASSO';
T = readtable(infile,'ReadVariableNames',false);
[y_lasso, idx_lasso] = ismember(T.Var1, sce.g);

infile = 'genes_all_qubo_manual_pseudo.txt';
highlightg = true;  ttxt = 'QUBO';
T = readtable(infile,'ReadVariableNames',false);
[y_qubo, idx_qubo] = ismember(T.Var1, sce.g);

idx_qubo = idx_qubo(1:20);
idx_lasso = idx_lasso(1:20);

assert( all(y_qubo) && all(y_lasso) )

X = sc_transform(sce.X);

g_lasso = sce.g(idx_lasso);
g_qubo = sce.g(idx_qubo);

X_lasso = X(idx_lasso,:);
X_qubo = X(idx_qubo,:);

[y,idx] = ismember('manual_pseudotime', sce.list_cell_attributes(1:2:end));
assert(all(y))
t = sce.list_cell_attributes{idx*2};

toc;

% Smoothing parameter
sp = 0.75;

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


%%

f = figure;
hold on
for k=1:length(g_qubo)
    h_lasso = plot(t_lasso_sort(:,k), y_lasso_fit(:,k),'LineWidth',2,'Color','k');
end
xlim([0 max(t)]);
ylim([-2 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('LASSO')
f.Position(3) = 545;
f.Position(4) = 266;

f=figure;
hold on
for k=1:length(g_qubo)
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',2,'Color', 'k');
end
xlim([0 max(t)]);
ylim([-3 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('QUBO')

[idx_common,~,c] = intersect(idx_lasso, idx_qubo);
for k=1:length(idx)
    if ismember(k,c)
        plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',0.1,'Color', [0, 0.8, 0]);
    end
end

f.Position(3) = 545;
f.Position(4) = 266;

fprintf('\nLASSO genes\n');
fprintf('%s, ', sort(sce.g(idx_lasso)));
fprintf('\nQUBO genes\n');
fprintf('%s, ', sort(sce.g(idx_qubo)));
fprintf('\nCommon genes\n');
fprintf('%s, ', sort(sce.g(idx_common)));
fprintf('\n');

%% Plotting with gene legends QUBO only

f = figure;
f.Position(3) = f.Position(3) * 1.4;
hold on;

colors = colormap("hsv");
ngq = length(g_qubo);

% Logical arrays for indices
idx_only_qubo = false(ngq, 1);
idx_common_qubo = false(ngq, 1);

% Compute indices without reordering
for ig = 1:ngq
    if ismember(g_qubo(ig), g_lasso)
        idx_common_qubo(ig) = true;
    else
        idx_only_qubo(ig) = true;
    end
end

% Plotting only unique qubo genes
legend_text = {};
h_lines = []; % To store plot handles

idx = find(idx_only_qubo); % Indices for unique genes
for i = 1:length(idx)
    gene_idx = idx(i);
    h_lines(i) = plot(t_qubo_sort(:, gene_idx), y_qubo_fit(:, gene_idx), ...
                      'LineWidth', 2, ...
                      'Color', colors(mod(gene_idx * 20 - 1, size(colors, 1)) + 1, :));
    legend_text{i} = g_qubo{gene_idx}; % Ensure proper cell array indexing
end

% Create legend
legend(h_lines, legend_text, 'Location', 'bestoutside');

xlim([0 max(t)]);
xlabel('Pseudotime');
ylabel('Standardized Expression');
box on;
title("Selected features - pseudotime prediction - QUBO only");
%% Plotting with gene legends LASSO only
f = figure;
f.Position(3) = f.Position(3) * 1.4;
hold on;

colors = colormap("hsv");
ngl = length(g_lasso);

% Logical arrays for indices
idx_only_lasso = false(ngl, 1);
idx_common_lasso = false(ngl, 1);

% Compute indices without reordering
for ig = 1:ngl
    if ismember(g_lasso(ig), g_qubo)
        idx_common_lasso(ig) = true;
    else
        idx_only_lasso(ig) = true;
    end
end

% Plotting only unique lasso genes
legend_text = {};
h_lines = []; % To store plot handles

idx = find(idx_only_lasso); % Indices for unique genes
for i = 1:length(idx)
    gene_idx = idx(i);
    h_lines(i) = plot(t_lasso_sort(:, gene_idx), y_lasso_fit(:, gene_idx), ...
                      'LineWidth', 2, ...
                      'Color', colors(mod(gene_idx * 20 - 1, size(colors, 1)) + 1, :));
    legend_text{i} = g_lasso{gene_idx}; % Ensure proper cell array indexing
end

% Create legend
legend(h_lines, legend_text, 'Location', 'bestoutside');

xlim([0 max(t)]);
xlabel('Pseudotime');
ylabel('Standardized Expression');
box on;
title("Selected features - pseudotime prediction - LASSO only");


%% Plotting with gene legends for overlapping genes (QUBO order)
f = figure;
f.Position(3) = f.Position(3) * 1.4;
hold on;

colors = colormap("hsv");
ngq = length(g_qubo);

% Logical array for indices in g_qubo
idx_common_qubo = false(ngq, 1);

% Compute indices for overlapping genes based on QUBO order
for ig = 1:ngq
    if ismember(g_qubo(ig), g_lasso)
        idx_common_qubo(ig) = true;
    end
end

% Plotting overlapping genes in QUBO order
legend_text = {};
h_lines = []; % To store plot handles

idx = find(idx_common_qubo); % Indices for overlapping genes in QUBO order
for i = 1:length(idx)
    gene_idx = idx(i);
    % Find corresponding index in LASSO
    lasso_idx = find(strcmp(g_lasso, g_qubo{gene_idx}));
    h_lines(i) = plot(t_lasso_sort(:, lasso_idx), y_lasso_fit(:, lasso_idx), ...
                      'LineWidth', 2, ...
                      'Color', colors(mod(gene_idx * 20 - 1, size(colors, 1)) + 1, :));
    legend_text{i} = g_qubo{gene_idx}; % Ensure proper cell array indexing
end

% Create legend
legend(h_lines, legend_text, 'Location', 'bestoutside');

xlim([0 max(t)]);
xlabel('Pseudotime');
ylabel('Standardized Expression');
box on;
title("Selected features - pseudotime prediction - Overlapping genes (QUBO order)");

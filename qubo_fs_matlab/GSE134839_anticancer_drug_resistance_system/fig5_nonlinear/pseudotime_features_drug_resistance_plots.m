load('../GSE134839_cleandata_lite.mat') % Original manuscript utilized this
g = sce.g;
X = full(sc_transform(sce.X, "type", "PearsonResiduals"));
K = 100;

%Tqubo = qfeatures_qubo_base(X, g, y, K, false);
g1 = Tqubo.selectedGenes;
%Tml = mlfeatures_base(X, g, y, K, 1);  
g2 = Tml{1}.selectedGenes;
% Random forest regression
g3 = Tml{4}.selectedGenes;

rfn_u_lasso = union(g2, g3);
qubo_only = setdiff(g1, rfn_u_lasso,'stable'); % This keeps g1 order
%qubo_only = g1;

qubo_shared = intersect(g1, g2, 'stable' );
writematrix(qubo_only','qubo_only_features.txt')
writematrix(qubo_shared','qubo_shared_features.txt')

idx_qubo = zeros(length(qubo_only), 1);
for i = 1:length(qubo_only)
    idx_qubo(i) = find(sce.g == qubo_only(i));
end
g_qubo = sce.g(idx_qubo);
X_qubo = X(idx_qubo,:);

idx_qubo_shared = zeros(length(qubo_shared), 1);
for i = 1:length(qubo_shared)
    idx_qubo_shared(i) = find(sce.g == qubo_shared(i));
end
g_qubo_shared = sce.g(idx_qubo_shared);
X_qubo_shared = X(idx_qubo_shared,:);

% Predictor
cell_type_target = 'manual_pseudotime';

% Preparing target predictor y from pseudo-time values per cell
idx = find(contains(sce.list_cell_attributes(1:2:end), cell_type_target));
if isempty(idx), return; end
t = sce.list_cell_attributes{idx*2};
toc;

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
% QUBO shared genes fitting 
ngene = size(X_qubo_shared, 1);
ncell = size(X_qubo_shared, 2);
y_qubo_shared_fit = zeros(ncell, ngene);
t_qubo_shared_sort = zeros(ncell, ngene);
for ig = 1:ngene
    [y_qubo_shared_fit(1:ncell, ig), idx] = ...
                        loess_smoothing(t, X_qubo_shared(ig,:)', sp);
    t_qubo_shared_sort(1:ncell, ig) = t(idx);
end
toc;

%%

f=figure;
hold on;
for k=1:length(g_qubo_shared)
    plot(t_qubo_shared_sort(:,k), y_qubo_shared_fit(:,k),'LineWidth',1,'Color', 'k');
end

hold on;
for k=1:length(qubo_only)
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',2,'Color', 'g');
end

xlim([0 max(t)]);
ylim([-2.5 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('QUBO')

%% Non-linear genes Possible figure 3
% Non-linear genes
%my_genes = g_qubo(1:10);
my_genes = g_qubo(21:28);

% Create a 2x5 figure
f=figure;

for idx = 1:length(my_genes)
    % subplot(2, 5, i); % Create subplot
    nexttile;
    hold on;

    % Plot all QUBO genes in gray
    for k = 1:length(qubo_only)
        plot(t_qubo_sort(:,k), y_qubo_fit(:,k), ...
            'LineWidth', 0.5, 'Color', 'k'); %[0.35, 0.35, 0.35]);
        %plot(t_qubo_sort(:,k), y_qubo_fit(:,k), 'LineWidth', 2, 'Color', 'k');

    end

    % Plot the selected gene in red
    gene_index = find(strcmp(g_qubo, my_genes(idx)));
    plot(t_qubo_sort(:,gene_index), ...
        y_qubo_fit(:,gene_index), 'LineWidth', 2, ...
        'Color', 'r');

    %xlim([-0.06 max(t)+0.04]);
    %ylim([-2 4]);
    
    %xlabel('Pseudotime')
    %ylabel('Standardized Expression')
    box on
    title(sprintf('%s', my_genes(idx)));
end
f.Position(3:4)=[1364 420];

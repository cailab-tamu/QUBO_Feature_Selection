tic;
load('../Data_hESC_EC_day1_5000g.mat')

infile = 'genes_all_rfr_splinefit.txt'; 
highlightg = false; ttxt = 'RFR';
T = readtable(infile,'ReadVariableNames',false);
[y_rfr, idx_rfr] = ismember(T.Var1, sce.g);

infile = 'genes_all_qubo_splinefit.txt';
highlightg = true;  ttxt = 'QUBO';
T = readtable(infile,'ReadVariableNames',false);
[y_qubo, idx_qubo] = ismember(T.Var1, sce.g);

assert( all(y_qubo) && all(y_rfr) )

X = sc_transform(sce.X);

g_rfr = sce.g(idx_rfr);
g_qubo = sce.g(idx_qubo);

X_rfr = X(idx_rfr,:);
X_qubo = X(idx_qubo,:);

[y,idx] = ismember('splinefit_pseudotime', sce.list_cell_attributes(1:2:end));
assert(all(y))
t = sce.list_cell_attributes{idx+1};

toc;

% Smoothing parameter
sp = 0.75;

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
for k=1:50
    h_rfr = plot(t_rfr_sort(:,k), y_rfr_fit(:,k),'LineWidth',2,'Color','k');
end
xlim([0 max(t)]);
ylim([-2 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('RFR')
f.Position(3) = 545;
f.Position(4) = 266;

f=figure;
hold on
for k=1:50
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',2,'Color', 'k');
end
xlim([0 max(t)]);
ylim([-2 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('QUBO')

[idx_common,~,c] = intersect(idx_rfr, idx_qubo);
for k=1:50
    if ismember(k,c)
        plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',0.1,'Color', [0, 0.8, 0]);
    end
end

f.Position(3) = 545;
f.Position(4) = 266;

fprintf('\nRFR genes\n');
fprintf('%s, ', sort(sce.g(idx_rfr)));
fprintf('\nQUBO genes\n');
fprintf('%s, ', sort(sce.g(idx_qubo)));
fprintf('\nCommon genes\n');
fprintf('%s, ', sort(sce.g(idx_common)));
fprintf('\n');

%% Non-linear genes Possible figure 3
my_genes = setdiff(g_qubo,g_rfr,'stable');
% Create a 2x5 figure
f=figure;

for idx = 1:length(my_genes)
    % subplot(2, 5, i); % Create subplot
    nexttile;
    hold on;

    % Plot all QUBO genes in gray
    for k = 1:50
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

%% Non-linear genes Possible figure 3
%my_genes = setdiff(g_qubo,g_rfr,'stable');
my_genes = setdiff(g_rfr,g_qubo,'stable');

% Create a 2x5 figure
f=figure;
for idx = 1:length(my_genes)
    % subplot(2, 5, i); % Create subplot
    nexttile;
    hold on;

    % Plot all QUBO genes in gray
    for k = 1:50
        plot(t_rfr_sort(:,k), y_rfr_fit(:,k), ...
            'LineWidth', 0.5, 'Color', 'k'); %[0.35, 0.35, 0.35]);
        %plot(t_qubo_sort(:,k), y_qubo_fit(:,k), 'LineWidth', 2, 'Color', 'k');

    end

    % Plot the selected gene in red
    gene_index = find(strcmp(g_rfr, my_genes(idx)));
    plot(t_rfr_sort(:,gene_index), ...
        y_rfr_fit(:,gene_index), 'LineWidth', 2, ...
        'Color', 'r');

    %xlim([-0.06 max(t)+0.04]);
    %ylim([-2 4]);
    
    %xlabel('Pseudotime')
    %ylabel('Standardized Expression')
    box on
    title(sprintf('%s', my_genes(idx)));
end
f.Position(3:4)=[1364 420];

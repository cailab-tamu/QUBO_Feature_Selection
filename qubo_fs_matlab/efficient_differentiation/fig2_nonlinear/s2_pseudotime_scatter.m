tic;
load('../Data_hESC_EC_day1_5000g.mat')

infile = 'genes_all_lasso_splinefit.txt'; 
highlightg = false; ttxt = 'LASSO';
T = readtable(infile,'ReadVariableNames',false);
[y_lasso, idx_lasso] = ismember(T.Var1, sce.g);

infile = 'genes_all_qubo_splinefit.txt';
highlightg = true;  ttxt = 'QUBO';
T = readtable(infile,'ReadVariableNames',false);
[y_qubo, idx_qubo] = ismember(T.Var1, sce.g);

assert( all(y_qubo) && all(y_lasso) )

X = sc_transform(sce.X);

g_lasso = sce.g(idx_lasso);
g_qubo = sce.g(idx_qubo);

X_lasso = X(idx_lasso,:);
X_qubo = X(idx_qubo,:);

[y,idx] = ismember('splinefit_pseudotime', sce.list_cell_attributes(1:2:end));
assert(all(y))
t = sce.list_cell_attributes{idx+1};

toc;

% Smoothin parameter
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

%% Plotting with Legend

%{
f = figure;
f.Position(3)=f.Position(3)*1.4;
hold on

my_color = [0, 0.8, 0];

h_lasso = plot(t_lasso_sort(:,1), y_lasso_fit(:,1),'LineWidth',3,'Color','k');
h_qubo = plot(t_qubo_sort(:,1), y_qubo_fit(:,1),'LineWidth',1,'Color', my_color);

% Plot remaining lines
for k = 2:size(y_lasso_fit,2)
    plot(t_lasso_sort(:,k), y_lasso_fit(:,k),'LineWidth',3,'Color','k');
end 
for k = 2:size(y_qubo_fit,2)
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',1,'Color', my_color);
end 

% Create legend
legend([h_lasso, h_qubo], {'LASSO', 'QUBO'},'Location','northwest');

xlim([0 max(t)]);
% ylim([-6 6]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title("Selected features - pseudotime prediction");
%}

%%

f = figure;
hold on
for k=1:35
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
for k=1:35
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',2,'Color', 'k');
end
xlim([0 max(t)]);
ylim([-2 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('QUBO')

[idx_common,~,c] = intersect(idx_lasso, idx_qubo);
for k=1:35
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

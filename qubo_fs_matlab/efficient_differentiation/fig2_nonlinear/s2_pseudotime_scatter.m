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
for k=1:50
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
for k=1:50
    plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',2,'Color', 'k');
end
xlim([0 max(t)]);
ylim([-2 4]);
xlabel('Pseudotime')
ylabel('Standardized Expression')
box on
title('QUBO')

[idx_common,~,c] = intersect(idx_lasso, idx_qubo);
for k=1:50
    if ismember(k,c)
        plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',0.1,'Color', [0, 0.8, 0]);
    end
end

% % Non-linear genes
% my_genes = ["IGFBP5" "KLK10" "MAP1B" "RGS10" "SFRP1" "TP53I11" "TRH" "TUBA1C" ...
%             "VIM" "YWHAB"]; 
% for k=1:50
%     if ismember(g_qubo(k), my_genes)
%         plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',0.1,'Color', [1, 0, 0]);
%     end
% end
f.Position(3) = 545;
f.Position(4) = 266;

fprintf('\nLASSO genes\n');
fprintf('%s, ', sort(sce.g(idx_lasso)));
fprintf('\nQUBO genes\n');
fprintf('%s, ', sort(sce.g(idx_qubo)));
fprintf('\nCommon genes\n');
fprintf('%s, ', sort(sce.g(idx_common)));
fprintf('\n');



%% Non-linear genes Possible figure 3
% Non-linear genes
my_genes = ["IGFBP5" "KLK10" "MAP1B" "RGS10" "SFRP1" "TP53I11" "TRH" "TUBA1C" "VIM" "YWHAB"]; 

% Create a 2x5 figure
figure;

for i = 1:length(my_genes)
    subplot(2, 5, i); % Create subplot
    hold on;

    % Plot all QUBO genes in gray
    for k = 1:50
        plot(t_qubo_sort(:,k), y_qubo_fit(:,k), 'LineWidth', 2, 'Color', [0.35, 0.35, 0.35]);
        %plot(t_qubo_sort(:,k), y_qubo_fit(:,k), 'LineWidth', 2, 'Color', 'k');

    end

    % Plot the selected gene in red
    gene_index = find(strcmp(g_qubo, my_genes(i)));
    plot(t_qubo_sort(:,gene_index), y_qubo_fit(:,gene_index), 'LineWidth', 1, 'Color', 'r');

    xlim([-0.06 max(t)+0.04]);
    ylim([-2 4]);
    %xlabel('Pseudotime')
    %ylabel('Standardized Expression')
    box on
    title(sprintf('%s', my_genes(i)));
end
%% Plotting with Legend per method

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

%% Plotting with gene legends QUBO only

% f = figure;
% f.Position(3)=f.Position(3)*1.4;
% hold on;
% 
% colors = colormap(jet);
% [g_common, idx_common_qubo] = intersect(g_qubo, g_lasso);
% [ g_qubo_unique, idx_only_qubo]= setdiff(g_qubo, g_common);
% 
% my_idx = [11, 12, 14, 19, 21, 26, 28, 29, 30];
% legend_text = {};
% h_lines = []; % To store plot handles
% 
% % Plot all lines and create legend entries
% ibeg = 11; 
% iend = 15;
% %for i = 1:length(idx_only_qubo)
% for i = ibeg:iend
%     k = idx_only_qubo(i);
%     h_lines(i) = plot(t_qubo_sort(:,k), y_qubo_fit(:,k), 'LineWidth', 1, 'Color', colors(mod(i*20-1, size(colors, 1)) + 1, :));
%     legend_text{i} = g_qubo_unique(i);
% end
% 
% % Create legend
% legend(h_lines, legend_text, 'Location', 'bestoutside');
% 
% xlim([0 max(t)]);
% xlabel('Pseudotime')
% ylabel('Standardized Expression')
% box on
% title("Selected features - pseudotime prediction");

%% Plotting with gene legend QUBO only in ranges
% 
% f = figure;
% f.Position(3)=f.Position(3)*1.4;
% hold on;
% colors = colormap(jet);
% [g_common, idx_common_qubo] = intersect(g_qubo, g_lasso);
% [ g_qubo_unique, idx_only_qubo]= setdiff(g_qubo, g_common);
% 
% legend_text = {};
% h_lines = []; % To store plot handles
% 
% ibeg = 11; 
% iend = 11;
% 
% % Plot all lines and create legend entries
% for i = ibeg:iend
%     k = idx_only_qubo(i);
%     h_lines(i-ibeg+1) = plot(t_qubo_sort(:,k), y_qubo_fit(:,k), 'LineWidth', 1, 'Color', colors(mod(i*40-1, size(colors, 1)) + 1, :));
%     legend_text{i-ibeg+1} = g_qubo_unique(i);
% end
% 
% % Create legend
% legend(h_lines, legend_text, 'Location', 'bestoutside');
% 
% xlim([0 max(t)]);
% xlabel('Pseudotime')
% ylabel('Standardized Expression')
% box on
% title("Selected features - pseudotime prediction");

rng default;
p = 10000;
n = 50;

source_f = [5,11,7,1,14];
target_f = [16,17,18,19,20];

% Source features are highly correlated to target features 
% where target is "regulated" source
[p,n,X,Y] = synthetic_data(p, n, source_f, target_f);

% Features to extract
K = 5;
g = 1:n;

X = X';
Y = Y';
Tqubo = qfeatures_qubo_base_v2(X, g, Y, K, false);

Tml  = mlfeatures_base_v2(X, g, Y, K, 1);

% % ------------------------- Plotting ------------------------------
% feat_lasso = Tml.sol_genes_lasso;
% X_lasso = X(feat_lasso,:);
% 
% feat_qubo = Tqubo.sol_genes;
% X_qubo = X(feat_qubo,:);
% 
% % Target predictor
% t = Y;
% 
% % Smoothing parameter
% sp = 0.75;
% 
% % LASSO fitting
% ngene = size(X_lasso, 1);
% ncell = size(X_lasso, 2);
% y_lasso_fit = zeros(ncell, ngene);
% t_lasso_sort = zeros(ncell, ngene);
% for ig = 1:ngene
%     [y_lasso_fit(1:ncell, ig), idx] = ...
%                         loess_smoothing(t, X_lasso(ig,:)', sp);
%     t_lasso_sort(1:ncell, ig) = t(idx);
% end
% 
% % QUBO fitting
% ngene = size(X_qubo, 1);
% ncell = size(X_qubo, 2);
% y_qubo_fit = zeros(ncell, ngene);
% t_qubo_sort = zeros(ncell, ngene);
% for ig = 1:ngene
%     [y_qubo_fit(1:ncell, ig), idx] = ...
%                         loess_smoothing(t, X_qubo(ig,:)', sp);
%     t_qubo_sort(1:ncell, ig) = t(idx);
% end
% 
% %% Plotting with Legend
% f = figure;
% f.Position(3)=f.Position(3)*1.4;
% hold on
% 
% my_color = [0, 0.8, 0];
% 
% h_lasso = plot(t_lasso_sort(:,1), y_lasso_fit(:,1),'LineWidth',3,'Color','k');
% h_qubo = plot(t_qubo_sort(:,1), y_qubo_fit(:,1),'LineWidth',1,'Color', my_color);
% 
% % Plot remaining lines
% for k = 2:size(y_lasso_fit,2)
%     plot(t_lasso_sort(:,k), y_lasso_fit(:,k),'LineWidth',3,'Color','k');
% end 
% for k = 2:size(y_qubo_fit,2)
%     plot(t_qubo_sort(:,k), y_qubo_fit(:,k),'LineWidth',1,'Color', my_color);
% end 
% 
% % Create legend
% legend([h_lasso, h_qubo], {'LASSO', 'QUBO'},'Location','northwest');
% 
% xlim([0 max(t)]);
% % ylim([-6 6]);
% xlabel('Predictor')
% ylabel('Standardized data')
% box on
% title("Selected features");
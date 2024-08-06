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
Tqubo = qfeatures_qubo_base(X, g, Y, K, false);

Tml  = mlfeatures_base(X, g, Y, K, 2);

inter_feat_qubo = intersect(Tqubo.sol_genes, source_f);
inter_feat_lasso = intersect(Tml.sol_genes_lasso, source_f);
inter_feat_relief = intersect(Tml.sol_genes_relief, source_f);

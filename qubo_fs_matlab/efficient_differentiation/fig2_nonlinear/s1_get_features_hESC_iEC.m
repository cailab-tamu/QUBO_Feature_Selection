load('../Data_hESC_EC_day1_5000g.mat') % Original manuscript utilized this
g = sce.g;
X = full(sc_transform(sce.X, "type", "PearsonResiduals"));
K = 50;

timetype = 'splinefit';

idx = find(contains(sce.list_cell_attributes(1:2:end), sprintf('%s_pseudotime',timetype)));
if isempty(idx), returen; end
y = sce.list_cell_attributes{idx*2};
y = y';

Tqubo = qfeatures_qubo_base(X, g, y, K, false);
g1 = Tqubo.sol_genes;

Tml = mlfeatures_base(X, g, y, K, 1);  
g2 = Tml.sol_genes_lasso;

qubo_n_lasso = intersect(g1, g2);
qubo_only = setdiff(g1, qubo_n_lasso);
lasso_only = setdiff(g2, qubo_n_lasso);

% run.web_Enrichr(g1)
% pause(5)
% run.web_Enrichr(g2)

%writetable(table([qubo_only; qubo_n_lasso]), sprintf('genes_picked_by_qubo_%s.txt', ...
%    timetype),'WriteVariableNames',false);

writetable(table(qubo_n_lasso'), sprintf('genes_shared_%s.txt', ...
    timetype),'WriteVariableNames',false);

writetable(table(qubo_only'), sprintf('genes_only_qubo_%s.txt', ...
    timetype),'WriteVariableNames',false);

writetable(table(lasso_only'), sprintf('genes_only_lasso_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');

writetable(table(Tqubo.sol_genes'), sprintf('genes_all_qubo_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');

writetable(table(Tml.sol_genes_lasso'), sprintf('genes_all_lasso_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');
    
numunique = length(qubo_only)



t1 = readtable('qubo_features.txt', 'ReadVariableNames', false);
t2 = readtable('lasso_features.txt', 'ReadVariableNames', false);
g1 = string(t1.Var1);
g2 = string(t2.Var1);

load('../GSE134839_cleandata_lite.mat') % Original manuscript utilized this
g = sce.g;
X = full(sc_transform(sce.X, "type", "PearsonResiduals"));

qubo_n_lasso = intersect(g1, g2);
qubo_only = setdiff(g1, qubo_n_lasso);
lasso_only = setdiff(g2, qubo_n_lasso);

% run.web_Enrichr(g1)
% pause(5)
% run.web_Enrichr(g2)

%writetable(table([qubo_only; qubo_n_lasso]), sprintf('genes_picked_by_qubo_%s.txt', ...
%    timetype),'WriteVariableNames',false);
timetype = 'manual_pseudo';

writetable(table(qubo_n_lasso'), sprintf('genes_shared_%s.txt', ...
    timetype),'WriteVariableNames',false);

writetable(table(qubo_only'), sprintf('genes_only_qubo_%s.txt', ...
    timetype),'WriteVariableNames',false);

writetable(table(lasso_only'), sprintf('genes_only_lasso_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');

writetable(table(g1), sprintf('genes_all_qubo_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');

writetable(table(g2), sprintf('genes_all_lasso_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');
    
numunique = length(qubo_only)


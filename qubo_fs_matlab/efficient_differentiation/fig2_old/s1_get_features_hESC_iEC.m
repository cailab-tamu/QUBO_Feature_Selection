load('../Data_hESC_EC_day1_5000g.mat') % Original manuscript utilized this
g = sce.g;
X = full(sc_transform(sce.X, "type", "PearsonResiduals"));
K = 50;

timetype = 'splinefit';

idx = find(contains(sce.list_cell_attributes(1:2:end), sprintf('%s_pseudotime',timetype)));
if isempty(idx), returen; end
y = sce.list_cell_attributes{idx*2};
y = y';

%Tqubo = qfeatures_qubo_base(X, g, y, K, false);
g1 = Tqubo.selectedGenes;

%Tml = mlfeatures_base(X, g, y, K, 'rfr');  
g2 = Tml{4}.selectedGenes;

qubo_n_rfr = intersect(g1, g2, 'stable');
qubo_only = setdiff(g1, qubo_n_rfr, 'stable');
rfr_only = setdiff(g2, qubo_n_rfr, 'stable');

% run.web_Enrichr(g1)
% pause(5)
% run.web_Enrichr(g2)

%writetable(table([qubo_only; qubo_n_rfr]), sprintf('genes_picked_by_qubo_%s.txt', ...
%    timetype),'WriteVariableNames',false);

writetable(table(qubo_n_rfr'), sprintf('genes_shared_%s.txt', ...
    timetype),'WriteVariableNames',false);

writetable(table(qubo_only'), sprintf('genes_only_qubo_%s.txt', ...
    timetype),'WriteVariableNames',false);

writetable(table(rfr_only'), sprintf('genes_only_rfr_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');

writetable(table(Tqubo.selectedGenes'), sprintf('genes_all_qubo_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');

writetable(table(Tml{4}.selectedGenes'), sprintf('genes_all_rfr_%s.txt', ...
    timetype),'WriteVariableNames',false,'Delimiter','\t');
    
numunique = length(qubo_only)


inter_qubo_rfr = intersect(Tqubo.selectedGenes,Tml{4}.selectedGenes)
g_nl = ["MAP1B" "KLK10" "TRH" "IGFBP5"  "TP53I11" "RGS10" "SFRP1" "TUBA1C" "VIM" "YWHAB"];

% compared to rfr
qubo_only_nl = setdiff(qubo_only, inter_qubo_rfr, 'stable');
qubo_only_nl = intersect(qubo_only_nl, g_nl,'stable');

% QUBO vs LASSO + RFR
g_lasso_n_rfr = union(Tml{4}.selectedGenes, Tml{1}.selectedGenes);
g_only_qubo = setdiff( Tqubo.selectedGenes, g_lasso_n_rfr,'stable');

writetable(table(g_only_qubo'), sprintf('genes_only_qubo_vs_all_%s.txt', ...
    timetype),'WriteVariableNames',false);
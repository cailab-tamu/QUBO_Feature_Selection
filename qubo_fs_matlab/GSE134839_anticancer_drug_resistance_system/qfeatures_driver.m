% Set path for Linux
%my_path  = "/scratch/user/ssromerogon/scGEAToolbox";
my_path = "../../src-v0.2";
addpath(genpath(my_path));

path="GSE134839_cleandata_lite.mat";
data  = load(path);
sce = data.sce;
clear data;

% Pre-processing 
ng = sce.NumGenes;
g = sce.g(1:ng);
X = full(sce.X(1:ng,:));
%X = sc_norm(X);
X = full(sc_transform(X, "type","PearsonResiduals"));

% Features to extract
K = 100; 

% Predictor
cell_type_target = "manual_pseudotime";

% Preparing target predictor y from pseudo-time values per cell
idx = find(contains(sce.list_cell_attributes(1:2:end), cell_type_target));
if isempty(idx), returen; end
y = sce.list_cell_attributes{idx*2};
%y = y';

fprintf("Final matrix size %d , %d \n",size(X));

% readR false will recompute R0 (MI)
readR = false;
Tqubo = qfeatures_qubo_base( X, g, y, K, readR);

modes = ["lasso", "elastic_net", "rrelieff", "fittree", "fsmrmr"];
Tml = cell(length(modes), 1);
for i = 1:length(modes)
    imode = modes(i); % Get the current mode
    fprintf("Current mode %s \n", imode); % Corrected variable name
    T = mlfeatures_base(X, g, y, K, imode); % Call the feature selection function
    T.mode = imode; % Add mode as a structure field
    Tml{i} = T; % Store the result
    outfile = sprintf("%s_ng_%d_features_%d.txt",imode,ng,K);
    writematrix(T.selectedGenes',outfile);
end

% Label stuff for saving tables
cts = strcat(string(K),"_",cell_type_target);
cts = strcat("_f",cts);
str_numcells=strcat("_cells_",string(sce.NumCells));
fname0 = strcat("_Genes_", int2str(ng));
fname0 = strcat(cts, fname0 );
fname0 = strcat(fname0, str_numcells); 

% Saving Tml
fname1 = strcat('Tml', fname0);
save(strcat(fname1,'.mat'),'Tml','-v7.3')

fname0 = strcat(fname0, str_numcells); 
fname1 = strcat('Tqubo', fname0);
save(strcat(fname1,'.mat'),'Tqubo','-v7.3')

writematrix(Tqubo.selectedGenes','qubo_features.txt');
writematrix(Tml{1}.selectedGenes',"lasso_features.txt");
writematrix(Tml{4}.selectedGenes',"fittree_features.txt");

% Energy landscape
load("R0.mat")
energy_path(R0, Tqubo.selectedGenes, Tml.selectedGenes, g, K,...
                 Tqubo.alphasol, 'energy_path.png');


% %% Saving matrices for d-wave
% load('R0.mat');
% 
% R = R0(1:end-1,1:end-1)/(K-1);
% J = R0(end,1:end-1);
% alpha = Tqubo.alphasol;
% 
% [~,qubo_sol ] = howmany(alpha,R,J);
% 
% % QUBO matrix (utilized for D-WAVE codes)
% Q = (1-alpha)*R - alpha*diag(J);
% writematrix(Q,'qubo_matrix.csv')
% writematrix(g,'genes.csv')

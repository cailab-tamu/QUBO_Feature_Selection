my_path = "../../qfeatures-src-v0.2_matlab";
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

modes = ["lasso", "elastic_net", "rrelieff", "fittree", "fsmrmr", "sequentialfs"];
Tml = cell(length(modes), 1);

for i = 1:length(modes)
    imoden = modes(i); % Get the current mode
    fprintf("Current mode %s \n", imoden); % Corrected variable name

    T = mlfeatures_base(X, g, Y, K, imoden); % Call the feature selection function

    T.mode = imoden; % Add mode as a structure field

    Tml{i} = T; % Store the result
end

% Label stuff for saving tables
fname0 = '_noisy_data'; 
% Saving Tml
fname1 = strcat('Tml', fname0);
save(strcat(fname1,'.mat'),'Tml','-v7.3')

fname1 = strcat('Tqubo', fname0);
save(strcat(fname1,'.mat'),'Tqubo','-v7.3')

inter_feat_qubo = intersect(Tqubo.selectedGenes, source_f)
inter_feat_lasso = intersect(Tml{1}.selectedGenes, source_f)
inter_feat_elastic_net = intersect(Tml{2}.selectedGenes, source_f)
inter_feat_relief = intersect(Tml{3}.selectedGenes, source_f)
inter_feat_fittree = intersect(Tml{4}.selectedGenes, source_f)
inter_feat_fsmrmr = intersect(Tml{5}.selectedGenes, source_f)
inter_feat_sequentialfs = intersect(Tml{6}.selectedGenes, source_f)




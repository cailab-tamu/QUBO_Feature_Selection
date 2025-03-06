my_path = "../../src-v0.2";
addpath(genpath(my_path));

path="GSE134839_cleandata_lite.mat";
data  = load(path);
sce = data.sce;
clear data;

% Pre-processing 
g = sce.g;
X = full(sce.X);
X = full(sc_transform(X, "type","PearsonResiduals"));

% Features to extract
K = 100; 

% Predictor
cell_type_target = "manual_pseudotime";

% Number of genes
ngenes = length(g);

% Preparing target predictor y from pseudo-time values per cell
idx = find(contains(sce.list_cell_attributes(1:2:end), cell_type_target));
if isempty(idx), returen; end
y = sce.list_cell_attributes{idx*2};
y = y';

fprintf("Final matrix size %d , %d \n",size(X));

[training_info, selectedGenes0] = cross_validation_qubo( X, g, y, K);
save('training_info.mat','training_info','-v7.3')
writematrix(selectedGenes0','qubo_features_train.txt');


load('training_info.mat')

% Initialize a cell array to store modified tables
modified_tables = cell(size(training_info));
for i = 1:numel(training_info)
    % Extract the current table
    tbl = training_info{i};
    
    % Remove 'Genes_trained' and 'Intersected_genes' columns
    tbl.Genes_trained = [];
    tbl.Intersected_genes = [];
    
    % Store the modified table
    modified_tables{i} = tbl;
end

% Concatenate the modified tables
combined_table = vertcat(modified_tables{:});

% Extract the data from the table
Itrain = combined_table{:,'I-train'};
TestAccuracyLocal = combined_table{:,'Test_accuracy_local(%)'};
TestAccuracy = combined_table{:,'Test_accuracy(%)'};

% Create a single figure for both plots
figure;

% Plot 1: Test_accuracy_local on the same figure
plot(Itrain, TestAccuracyLocal, 'o-', 'LineWidth', 1.5, 'Color', [0.2, 0.6, 0.2]); % Green color
hold on; % Keep the same plot

% Plot 2: Test_accuracy on the same figure
plot(Itrain, TestAccuracy, 's-', 'LineWidth', 1.5, 'Color', [0.2, 0.2, 0.8]); % Blue color

% Labels and title
xlabel('I-train');
ylabel('Accuracy (%)');
title(sprintf('I-train vs Test Accuracy for %d Features', K));

% Add a legend
%legend('Test Accuracy Local', 'Test Accuracy', 'Location', 'best');
legend('Test Accuracy Local', 'Test Accuracy', 'Location', 'southeast');

% Add grid for better readability
grid on;

% Hold off to stop adding plots to the same figure
hold off;

saveas(gcf, 'cross_validation_k100', 'svg');  % Save the figure to the specified path

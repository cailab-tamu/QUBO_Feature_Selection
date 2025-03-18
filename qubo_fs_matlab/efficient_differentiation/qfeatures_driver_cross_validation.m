my_path = "../../qfeatures-src-v0.2_matlab";
addpath(genpath(my_path));

path="Data_hESC_EC_day1_5000g.mat";
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
cell_type_target = "monocle3_pseudotime";

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
f = figure;

% Plot 1: Test_accuracy_local on the same figure
plot(Itrain, TestAccuracyLocal, 'o-', 'LineWidth', 1.5, 'Color', [0.2, 0.6, 0.2]); % Green color
hold on; % Keep the same plot

% Plot 2: Test_accuracy on the same figure
plot(Itrain, TestAccuracy, 's-', 'LineWidth', 1.5, 'Color', [0.2, 0.2, 0.8]); % Blue color

% Labels and title
xlabel('I-train','FontSize',14);
ylabel('Accuracy (%)','FontSize',14);
title(sprintf('I-train vs test accuracy for %d features', K),'FontSize',16);
ylim([90 100])

% Add a legend
lgd = legend('Test accuracy local', 'Test accuracy', 'Location', 'best');
%lgd = legend('Test Accuracy Local', 'Test Accuracy', 'Location', 'southeast');
fontsize(lgd, 11,'points')

% Add grid for better readability
grid on;

% Hold off to stop adding plots to the same figure
hold off;


% Save as high-quality PNG using print
filename = 'cross_validation_k100.png';
resolution = 300; 
print(f, filename, '-dpng', sprintf('-r%d', resolution));

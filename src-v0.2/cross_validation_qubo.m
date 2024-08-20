function [training_info, selectedGenes0, avg_training_accu] = cross_validation_qubo(X, g, y, K)
                    
    nFolds = 10; % Number of cross-validation folds
    cv = cvpartition(size(y, 2), 'KFold', nFolds);
    disp(cv)
    readR = false; % We need to test several times the MI for test and train sets

    fprintf('Running full set...\n');
    % Run QUBO feature selection on the full dataset
    [Tqubo_train, ~] = qfeatures_qubo_base(X, g, y, K, readR);
    selectedGenes0 = Tqubo_train.sol_genes;
    
    % Store selected genes for each fold
    training_info = cell(nFolds, 1);

    avg_training_accu = 0.0;
    for i = 1:nFolds
        fprintf('Running fold %d/%d\n', i, nFolds);
        
        % Partition data into training and test sets for this fold
        trainIdx = cv.training(i);
        testIdx = cv.test(i);

        X_train = X(:, trainIdx);
        y_train = y(trainIdx);
        
        X_test = X(:, testIdx);
        y_test = y(testIdx);

        % Run QUBO feature selection on the training set
        fprintf('X_train size ( %d, %d) \n ', size(X_train, 1), size(X_train, 2));
        [Tqubo_train, sa_train] = qfeatures_qubo_base(X_train, g, y_train, K, readR);

        % Store the selected genes for this fold
        genesTraining = Tqubo_train.sol_genes;

        % Evaluate the model on the test set
        [~, test_accu] = ifold_test(Tqubo_train, sa_train, X_test, y_test, K);
        fprintf("*****Test Accuracy (%%): %f \n", test_accu);

        avg_training_accu = test_accu + avg_training_accu;
        % Compute intersection of each fold's selected genes with selectedGenes0
        inter_genes = intersect(genesTraining, selectedGenes0);
        ninter = length(inter_genes);
        if isempty(inter_genes)
            fprintf('No intersected genes in fold %d\n', i);
        else
            fprintf('Number of intersected genes in fold %d: %d\n', i, ninter);
        end

        % Store fold-specific information
        training_info{i} = table(i, genesTraining, Tqubo_train.alphasol, ...
                                 inter_genes, ninter, test_accu);
        training_info{i}.Properties.VariableNames = {'I-train', 'Genes_trained',...
                                                     'Alpha_trained', 'Intersected_genes', ...
                                                     'No_intersected_genes', 'Test_accuracy(%)'};
    end
    avg_training_accu = avg_training_accu / nFolds;
    fprintf("Average training accuracy : %f (%%) \n", avg_training_accu);
end
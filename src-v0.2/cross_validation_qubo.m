function [errors, selectedGenesList, selectedGenes0, ...
          eners, ener_per_fit_all, intersections] = cross_validation_qubo(X, g, y, K)
                    
    nFolds = 10; % Number of cross-validation folds
    cv = cvpartition(size(y, 2), 'KFold', nFolds);
    disp(cv)
    selectedGenesList = cell(nFolds, 1); % Store selected genes for each fold
    intersections = cell(nFolds, 1); % Store intersections for each fold
    errors = zeros(nFolds + 1, 1);
    eners = zeros(nFolds + 1, 1);
    ener_per_fit_all = zeros(nFolds + 1, K);
    for i = 1:nFolds
        fprintf('Running fold %d/%d\n', i, nFolds);
        trainIdx = cv.training(i);
        testIdx = cv.test(i);
    
        X_train = X(:, trainIdx);
        y_train = y(trainIdx);
    
        X_test = X(:, testIdx);
        y_test = y(testIdx);
    
        % Run QUBO feature selection on the training set
        readR = false;
        fprintf('X_train size ( %d, %d) \n ', size(X_train, 1), size(X_train, 2));
        [Tqubo, sa_sol] = qfeatures_qubo_base(X_train, g, y_train, K, readR);

        % Store the selected genes for this fold
        selectedGenesList{i} = Tqubo.sol_genes;
    
        % Obtain QUBO function value, build energy landscape
        [eners(i), ener_per_fit_all(i,:), errors(i)] =  ...
            ifold_test(Tqubo, sa_sol, X_test, y_test, g, K);

    end

    fprintf('Running full set...\n');
    % Run QUBO feature selection on the full dataset
    [Tqubo, sa_sol] = qfeatures_qubo_base(X, g, y, K, readR);
    selectedGenes0 = Tqubo.sol_genes;

    % Obtain QUBO function value, build energy landscape
    [eners(nFolds+1), ener_per_fit_all(nFolds+1,:), errors(nFolds+1)] = ...
                        ifold_test(Tqubo, sa_sol, X, y, g, K);

    % Compute intersection of each fold's selected genes with selectedGenes0
    for i = 1:nFolds
        intersections{i} = intersect(selectedGenesList{i}, selectedGenes0);
        if isempty(intersections{i})
            fprintf('No intersected genes in fold %d\n', i);
        else
            fprintf('Number of intersected genes in fold %d: %d\n', i, length(intersections{i}));
        end
    end
end
function [Tsol] = mlfeatures_base(X, g, y, K, imode, alpha)
    % mlfeatures_base computes the feature selection (FS)
    % from count matrix X, genes g and y target.
    % INPUT:
    % X =====> Single cell count matrix
    % g =====> Genes/features from single cell experiment
    % y =====> Predictor/target variable
    % K =====> Number of features to retrieve
    % imode ==> imode: 1 compute LASSO FS
    %          imode: 2 compute Elastic Net FS
    %          imode: 3 compute ReliefF FS
    %          imode: 4 compute fitrtree FS
    %          imode: 5 compute sequentialfs FS
    %          imode: 6 compute fsrmrmr FS
    % alpha ==> Elastic net parameter
    % htos ==> Highest Time of Selection (for sequentialfs)
    % OUTPUT: 
    % Tsol ==> MATLAB table containing features and computation time

    % imode can be optional
    if nargin < 5; imode = "lasso"; end
    if nargin < 6; alpha = 0.5; end

    % Ensure K does not exceed the number of features
    if K > size(X, 1)
        error('K cannot be larger than the number of features.');
    end

    % Initialize empty variables to avoid undefined variable errors
    sol_genes = {};
    abs_coef = []; % Initialize abs_coef

    % Convert from sparse to full, if necessary
    if issparse(X)
        X = full(X);
    end
    if issparse(y)
        y = full(y);
    end

    % Start feature selection
    tic;
    options = statset('UseParallel',true);
    switch imode
        case "lasso"
            disp('LASSO feature selection activated');
            [B, FitInfo] = lasso(X', y, 'CV', 10, 'Options', options);
            idxLambda1SE = FitInfo.Index1SE;
            coef = B(:, idxLambda1SE);
            [abs_coef, sorted_idx] = sort(abs(coef), 'descend');
            abs_coef = abs_coef(1:K);
            selectedFeatures = sorted_idx(1:K);
            sol_genes = g(selectedFeatures);

        case "elastic_net"
            disp('Elastic Net feature selection activated');
            [B, FitInfo] = lasso(X', y, 'CV', 10, 'Alpha', alpha,'Options', options);
            idxLambda1SE = FitInfo.Index1SE;
            coef = B(:, idxLambda1SE);
            [abs_coef, sorted_idx] = sort(abs(coef), 'descend');
            abs_coef = abs_coef(1:K);
            selectedFeatures = sorted_idx(1:K);
            sol_genes = g(selectedFeatures);

        case "rrelieff"
            disp('ReliefF feature selection activated');
            nearest_n = 10;
            [selectedFeatures, coef] = relieff(X', y', nearest_n);
            selectedFeatures = selectedFeatures(1:K);
            abs_coef = abs(coef(selectedFeatures));
            sol_genes = g(selectedFeatures);

        case "fittree"
            disp('fitrensemble (Tree) feature selection activated');
            Mdl = fitrensemble(X', y, 'Method', 'Bag', 'NumLearningCycles', 100, 'Options', options);
            featureImportance = oobPermutedPredictorImportance(Mdl,'Options',options);
            abs_featureImportance = abs(featureImportance); % Take the absolute value
            [abs_coef, sorted_idx] = sort(abs_featureImportance, 'descend');
            abs_coef = abs_coef(1:K);
            selectedFeatures = sorted_idx(1:K);
            sol_genes = g(selectedFeatures);

        case "fsmrmr"
            disp('fsrmrmr feature selection activated');
            [selectedFeatures, coef] = fsrmrmr(X', y'); % Remove 'NumFeatures'
            selectedFeatures = selectedFeatures(1:K); % Select top K features
            abs_coef = abs(coef(selectedFeatures));
            sol_genes = g(selectedFeatures);

        case "sequentialfs"
            disp('sequentialfs feature selection activated');
            fun = @(Xtrain, ytrain, Xtest, ytest) loss(fitrtree(Xtrain, ytrain), Xtest, ytest);
            [selectedFeatures, ~] = sequentialfs(fun, X', y', 'CV', 10, 'nfeatures', K, ...
                                              'Options', options, 'direction', ...
                                              'forward');
            selectedFeatures = find(selectedFeatures==1);
            sol_genes = g(selectedFeatures);
            abs_coef = ones(1,length(selectedFeatures)); % sequentialfs doesnt return feature importance

    end
    time = toc;
    fprintf("FS time: %f \n", time);

    % Check the non-zero features (for cases where it makes sense)
    if ~isempty(abs_coef)
        idx = abs_coef > 0;
        sol_genes = sol_genes(idx);
        selectedFeatures = selectedFeatures(idx);
        abs_coef = abs_coef(idx);
    else
        idx = true(size(sol_genes)); % For sequentialfs or fsrmrmr, keep all
    end

    nsel = sum(idx);
    if nsel < K && imode ~="fsmrmr" && imode ~= "sequentialfs"  %dont warn on sequentialfs or fsrmrmr
        fprintf("Not able to select %d genes, providing %d \n", K, nsel);
    end

    % Ensure sol_genes is in the correct orientation
    if size(sol_genes, 1) > 1
        sol_genes = sol_genes';
    end
    if size(selectedFeatures, 1) > 1
        selectedFeatures = selectedFeatures';
    end
    if size(abs_coef, 1) > 1
        abs_coef = abs_coef';
    end

    % Create output table
    Tsol = table(sol_genes, selectedFeatures, abs_coef, time, ...
        'VariableNames', {'selectedGenes', 'featureIndices', ...
        'abs_coef', 'computationTime'});
end
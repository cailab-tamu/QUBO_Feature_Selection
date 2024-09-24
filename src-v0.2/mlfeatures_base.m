function [Tsol] = mlfeatures_base(X, g, y, K, mode, alpha)
    % mlfeatures_base computes the feature selection (FS)
    % from count matrix X, genes g and y target.
    % INPUT:
    % X =====> Single cell count matrix
    % g =====> Genes/features from single cell experiment
    % y =====> Predictor/target variable
    % K =====> Number of features to retrieve
    % mode ==> mode: 1 compute LASSO FS
    %          mode: 2 compute Elastic Net FS
    %          mode: 3 compute ReliefF FS
    % OUTPUT: 
    % Tsol ==> MATLAB table containing features and computation time

    % mode can be optional
    if nargin < 5; mode = 1; end
    if nargin < 6; alpha = 0.5; end

    % Ensure K does not exceed the number of features
    if K > size(X, 1)
        error('K cannot be larger than the number of features.');
    end

    % Initialize empty variables to avoid undefined variable errors
    sol_genes = {};

    % Convert from sparse to full, if necessary
    if issparse(X)
        X = full(X);
    end
    if issparse(y)
        y = full(y);
    end

    % Start feature selection
    tic;
    switch mode 
        case 1
            disp('LASSO feature selection activated');
            % LASSO solver with cross validation (CV) 10
            [B, FitInfo] = lasso(X', y, 'CV', 10);
            % Select features within 1 standard error (1SE)
            idxLambda1SE = FitInfo.Index1SE;
            % Score value from LASSO
            coef = B(:, idxLambda1SE);
            % Select top K features
            [~, sorted_idx] = sort(abs(coef), 'descend');
            selectedFeatures = sorted_idx(1:K);
            sol_genes = g(selectedFeatures);

        case 2
            disp('Elastic Net feature selection activated');
            % Elastic Net solver with cross validation (CV) 10
            [B, FitInfo] = lasso(X', y, 'CV', 10, 'Alpha', alpha);
            idxLambda1SE = FitInfo.Index1SE;
            coef = B(:, idxLambda1SE);
            [~, sorted_idx] = sort(abs(coef), 'descend');
            selectedFeatures = sorted_idx(1:K);
            sol_genes = g(selectedFeatures);

        case 3 
            disp('ReliefF feature selection activated');
            nearest_n = 10;
            [rankedFeatures, ~] = relieff(X', y', nearest_n);
            selectedFeatures = rankedFeatures(1:K);
            sol_genes = g(selectedFeatures);
    end
    time = toc;
    fprintf("FS time: %f \n", time);

    % Ensure sol_genes is in the correct orientation
    if size(sol_genes, 1) > 1
        sol_genes = sol_genes';
    end
    if size(selectedFeatures, 1) > 1
        selectedFeatures = selectedFeatures';
    end

    % Create output table
    Tsol = table(sol_genes, selectedFeatures, time, ...
         'VariableNames', {'selectedGenes', 'featureIndices', 'computationTime'});
end


function [Tsol] = mlfeatures_base(X, g, y, K, mode)
    % mlfeatures_base computes the feature selection (FS)
    % from count matrix X, genes g and y target.
    % INPUT:
    % X =====> Single cell coung matrix
    % g =====> Genes/features from single cell experiment
    % y =====> Predictor/target variable
    % K =====> Number of features to retrieve
    % mode ==> mode: 1 compute LASSO FS
    %          mode: 2 compute LASSO and relieff FS
    % OUTPUT: 
    % Tsol ==> MATLAB table containing features and computation time

    % mode can be optional
    if nargin < 5
        mode = 1;
    end

    % Initialize empty variables to avoid undefined variable errors
    sol_genes_lasso = {};
    sol_genes_relief = {};
    time_lasso = NaN;
    time_relief = NaN;

    % Conver from sparse to full
    X = full(X);
    y = full(y);

    % Lasso feature selection
    if mode >= 1
        tic;
        % LASSO solver with cross validation (CV) 10
        [B, FitInfo] = lasso(X', y, 'CV', 10);
        % Select features within 1 standard error (1SE)
        idxLambda1SE = FitInfo.Index1SE;
        % Score value from LASSO
        coef = B(:, idxLambda1SE);
        % Select top K features
        [~, sorted_idx] = sort(abs(coef), 'descend');
        selectedFeatures = sorted_idx(1:K);
        sol_genes_lasso = g(selectedFeatures);
        nx = size(sol_genes_lasso, 1);
        if nx > 1
            sol_genes_lasso = sol_genes_lasso';
        end
        time_lasso = toc;
        fprintf("LASSO FS time: %f \n", time_lasso);
    end
    
    % ReliefF feature selection
    if mode >= 2
        tic;
        % 10 nearest neighboars into account 
        nearest_n = 10;
        % relieff solver
        [rankedFeatures, ~] = relieff(X', y', nearest_n);
        % Select top K features
        sol_genes_relief = g(rankedFeatures(1:K));
        nx = size(sol_genes_relief, 1);
        if nx > 1
            sol_genes_relief = sol_genes_relief';
        end
        time_relief = toc;
        fprintf("ReliefF FS time: %f \n", time_relief);
    end

    % Create output table based on the selected mode
    switch mode
        case 1
            Tsol = table(sol_genes_lasso, time_lasso);
        case 2
            Tsol = table(sol_genes_lasso, sol_genes_relief, ...
                         time_lasso, time_relief);
        otherwise
            error('Invalid mode. Mode should be 1 or 2.');
    end
end

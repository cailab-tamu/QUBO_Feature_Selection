function [eners, ener_per_fit_all, errors] = ifold_test(Tqubo, sa_sol, ...
                                               X, y_test, g, K)
    % Obtain QUBO function value, build energy landscape
    % Loading MI
    load('R0.mat')
    % Re-built R (redundancy) and J (importance)
    R = R0(1:end-1, 1:end-1) / (K - 1);
    J = R0(end, 1:end-1);
    % Re-compute Q  (balance of R and J)
    Q = (1 - Tqubo.alphasol) * R - Tqubo.alphasol * diag(J);
    
    % Function value per feature accoring best solution and Q matrix
    ener_per_feat = Q * sa_sol.BestX;
    % Obtain top most important features
    [ener_per_feat, sort_idx] = sort(ener_per_feat, 'ascend');
    jdx = sa_sol.BestX(sort_idx) == 1;
    ener_per_fit_all = ener_per_feat(jdx);
    eners = Tqubo.fval;
    
    % Evaluate the model on the full dataset
    X_test_selected = X(ismember(g, Tqubo.sol_genes), :)';
    mdl = fitlm(X_test_selected, y_test);
    errors = mdl.RMSE;
    % mdl = fitlm(X_test_selected, y); % or this
    fprintf('Error %f \n', errors);

end
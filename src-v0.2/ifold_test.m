function [ener_test, training_accu] = ifold_test(Tqubo, sa_sol, X_test, y_test, K)
    % Obtain QUBO function value, build energy landscape
    % Loading MI

    tic;
    % Construct MI for test set
    data = [X_test; y_test];
    R0 = MI_construction(data);
    time_mi = toc;
    fprintf("MI construction (ifold_test) time: %f \n", time_mi);

    %fprintf("Looking for %d genes \n", K);
    % Redundancy matrix
    R = R0(1:end-1, 1:end-1) / (K - 1);
    % Importance vector
    J = R0(end, 1:end-1);

    % Re-compute Q  (balance of R and J) with trained alpha
    Q = (1 - Tqubo.alphasol) * R - Tqubo.alphasol * diag(J);
    
    % Function value per feature accoring training solution and Q matrix
    ener_per_feat = Q * sa_sol.BestX;

    % Obtain top most important features
    % [ener_per_feat, sort_idx] = sort(ener_per_feat, 'ascend');
    % jdx = sa_sol.BestX(sort_idx) == 1;
    % ener_per_g_test = ener_per_feat(jdx);

    % Energy value of training solution in test set
    ener_test = sa_sol.BestX' *ener_per_feat;
    % Error within test set and training set
    error_abs_pct = abs( (ener_test - Tqubo.fval) / Tqubo.fval);
    training_accu = 100*(1 - error_abs_pct);  

    % mdl = fitlm(X_test_selected, y); % or this
    %fprintf('Training Accuracy (%%) %f  \n', training_accu);

end
function [ener_accu, ener_embed_accu] = ifold_test(Tqubo_train, sa_train, ener0, Q);
    % ifold_test evaluates the trained solution and evaluates it in the 
    % test set and in the whole problem Q.
    % INPUTS:
    % Tqubo_train =====> training info such as alpha and cost function value
    % sa_train ========> solution vector from training
    % sa0 =============> solution of the whole problem (Q)
    % Q ===============> the full problem
    % OUTPUT:
    % ener_accu =======> energy evaluate into test Q

    % Test Q' solution x' in the original cost function Q
    evec = Q * sa_train.BestX;
    ener_train = sa_train.BestX' * evec;

    % Local energy solution (Q') for the training set
    ener_local = Tqubo_train.fval;

    % Measure error of trained solution in original Q compared to full solution
    error_abs_pct = abs( (ener_train - ener0) / ener0);
    ener_embed_accu = 100*(1 - error_abs_pct);  

    % Measure error of test solution compared to full solution
    error_abs_pct = abs( (ener_local - ener0) / ener0);
    ener_accu = 100*(1 - error_abs_pct);  

end
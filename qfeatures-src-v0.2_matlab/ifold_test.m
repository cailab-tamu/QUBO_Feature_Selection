function [ener_accu, ener_embed_accu] = ifold_test(Tqubo_train, sa_train, sa0, Q);
    % Obtain QUBO function value, build energy landscape
    % Loading MI

    % Asses the CV selection in the original cost function Q
    evec = Q * sa_train.BestX;
    ener_train = sa_train.BestX' * evec;
 
    % Asses the problem (all obs) selection in the original cost function Q
    evec = Q * sa0.BestX;
    ener0 = sa0.BestX' * evec;

    % Local energy solution for the training set
    ener_local = Tqubo_train.fval;

    % Measure deviation from annealed solution embedded in the original
    % cost function 
    error_abs_pct = abs( (ener_train - ener0) / ener0);
    ener_embed_accu = 100*(1 - error_abs_pct);  

    % Measure deviation from annealed solution embedded in the original
    % cost function 
    error_abs_pct = abs( (ener_local - ener0) / ener0);
    ener_accu = 100*(1 - error_abs_pct);  

end
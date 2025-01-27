function energy_path(R0, qubo_genes, lasso_genes, genes, K, alphasol, save_path)
    % Calculate the energy landscape for QUBO and LASSO solutions
    % Inputs:
    % R0: Initial redundancy matrix with importance vector as the last row
    % qubo_genes: Selected genes from QUBO solution
    % lasso_genes: Selected genes from LASSO solution
    % genes: Complete list of genes
    % K: Number of genes
    % alphasol: Trained alpha solution
    % save_path: Path to save the image (optional)

    % Redundancy matrix (scaled)
    R = R0(1:end-1, 1:end-1) / (K - 1);

    % Importance vector
    J = R0(end, 1:end-1);

    % Recompute Q matrix (balance of R and J) using trained alpha
    Q = (1 - alphasol) * R - alphasol * diag(J);

    % QUBO solution vector and energy calculation
    sol_qubo = ismember(genes, qubo_genes);
    ener_per_feat_qubo = Q * sol_qubo; % Compute energy for each feature
    ener_per_feat_qubo = sort(ener_per_feat_qubo); % Sort energies

    % LASSO solution vector and energy calculation
    sol_lasso = ismember(genes, lasso_genes);
    ener_per_feat_lasso = Q * sol_lasso; % Compute energy for each feature
    ener_per_feat_lasso = sort(ener_per_feat_lasso); % Sort energies

    % Cumulative energy calculations
    nsolq = min(K,length(qubo_genes));
    sum_ener_qubo = cumsum(ener_per_feat_qubo(1:nsolq));
    nsoll = min(K,length(lasso_genes));
    sum_ener_lasso = cumsum(ener_per_feat_lasso(1:nsoll));
    
    % Plot the energy landscapes
    range1 = 1:nsolq;

    % Plot lines with colors
    plot(range1, sum_ener_qubo, '-g', 'DisplayName', 'QUBO'); % Green line
    hold on;
    range2 = 1:nsoll;
    plot(range2, sum_ener_lasso, '-k', 'DisplayName', 'LASSO'); % Black line

    % Add circles to each data point
    plot(range1, sum_ener_qubo, 'og', 'MarkerSize', 3,'MarkerFaceColor', 'g', ...
         'DisplayName', 'QUBO Points'); % Green circles
    plot(range2, sum_ener_lasso, 'ok','MarkerSize', 3, 'MarkerFaceColor', 'k', ...
        'DisplayName', 'LASSO Points'); % Black circles 

    % Labels and legend
    xlabel('Number of Features');
    ylabel('Cumulative Energy');
    legend show;
    title('Energy Path Comparison');
    grid on;
    hold off;

    % Save the plot to a file if save_path is provided
    if nargin > 6 && ~isempty(save_path)
        saveas(gcf, save_path, 'svg');  % Save the figure to the specified path
        fprintf('Plot saved to %s\n', save_path);
    end
end

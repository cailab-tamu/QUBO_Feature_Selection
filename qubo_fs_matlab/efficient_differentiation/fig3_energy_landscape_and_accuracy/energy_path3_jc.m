function energy_path3_jc(R0, qubo_genes, lasso_genes, rfr_genes, genes, K, alphasol, save_path)
    % Calculate the energy landscape for QUBO, LASSO, and RFR solutions
    % Inputs:
    % R0: Initial redundancy matrix with importance vector as the last row
    % qubo_genes: Selected genes from QUBO solution
    % lasso_genes: Selected genes from LASSO solution
    % rfr_genes: Selected genes from RFR solution
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
    % RFR solution vector and energy calculation
    sol_rfr = ismember(genes, rfr_genes);
    ener_per_feat_rfr = Q * sol_rfr; % Compute energy for each feature
    ener_per_feat_rfr = sort(ener_per_feat_rfr); % Sort energies
    % Cumulative energy calculations
    nsolq = min(K,length(qubo_genes));
    sum_ener_qubo = cumsum(ener_per_feat_qubo(1:nsolq));
    nsoll = min(K,length(lasso_genes));
    sum_ener_lasso = cumsum(ener_per_feat_lasso(1:nsoll));
    nsolr = min(K,length(rfr_genes));
    sum_ener_rfr = cumsum(ener_per_feat_rfr(1:nsolr));
    
    % Plot the energy landscapes
    range1 = 1:nsolq;
    f = figure;
    % Plot lines with colors
    
    
    
    
    range3 = 1:nsolr;
    plot(range3, sum_ener_rfr, '-b', 'LineWidth', 2, 'DisplayName', 'RFR'); % Blue line
    hold on;
    range2 = 1:nsoll;
    plot(range2, sum_ener_lasso, '-r', 'LineWidth', 2, 'DisplayName', 'LASSO'); % Black line
    plot(range1, sum_ener_qubo, '-g', 'LineWidth', 2, 'DisplayName', 'QUBO'); % Green line

    xlim([0 330]);
    % Labels and legend
    xlabel('Number of Features','FontSize', 14);
    ylabel('Cumulative Energy','FontSize', 14);
    lgd = legend('show', 'Location','southwest');
    fontsize(lgd, 11,'points')
    title('Energy Path','FontSize', 15);
    grid on;
    hold off;
    % Save the plot to a file if save_path is provided
    if nargin > 7 && ~isempty(save_path)
        resolution = 300; 
        print(f, save_path, '-dpng', sprintf('-r%d', resolution));
        fprintf('Plot saved to %s\n', save_path);
    end
end
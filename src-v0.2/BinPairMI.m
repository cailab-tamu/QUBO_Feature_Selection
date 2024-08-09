function MI = BinPairMI(x, y)
    % BinPairMI computes the mutual information within x and y
    % INPUT:
    % x ====> Feature x observations
    % y ====> Feature y observations
    % OUTPUT
    % MI ===> Mutual information value within x-y

    % Estimate joint probability distribution with MATLAB automatic binning
    [joint_counts, ~, ~] = histcounts2(x, y);
    joint_prob = joint_counts ./ sum(joint_counts(:));
    
    % Estimate marginal probabilities
    x_marginal = sum(joint_prob, 2);
    y_marginal = sum(joint_prob, 1);
    
    % Handle zeros (add eps for numerical stability)
    eps0 = eps(realmin('single'));
    joint_prob = joint_prob + eps0;
    x_marginal = x_marginal + eps0;
    y_marginal = y_marginal + eps0;
    
    % Calculate entropy terms
    entropy_xy = -sum(joint_prob(:) .* log2(joint_prob(:)) );
    entropy_x = -sum(x_marginal .* log2(x_marginal) );
    entropy_y = -sum(y_marginal .* log2(y_marginal) );
    
    % Mutual information
    MI = entropy_x + entropy_y - entropy_xy;
end
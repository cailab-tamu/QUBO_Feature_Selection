function [p,n,X,y] = synthetic_data(p, n, sorce, target)
    nsource = size(sorce, 2);
    B = randn(n,n);
    R = corrcov(B'*B);
    X = mvnrnd(zeros(n,1),R,p);
    % Source features are highly correlated to target features 
    corrStd = 0.1;
    X(:,target) = X(:,sorce) + corrStd*randn(p,nsource);
    noiseStd = 0.1;
    % % Original function
    % t = 0.5*cos( X(:,sorce(4)) ) + ...
    %     sin( X(:,sorce(3)).*X(:,sorce(2)) ) + ...
    %     0.1*exp(X(:,sorce(5))).*log2(abs(X(:,sorce(1)))) + ...
    %     noiseStd*randn(p,1);
    %Original function
    t = 0.5*cos( 7*X(:,sorce(4)) ) + ... % 7 to denoise the randomness
        sin( X(:,sorce(3)).*X(:,sorce(2)) ) + ...
        0.1*exp(X(:,sorce(5))).*log2( abs(10*X(:,sorce(1))) ) + ...
        noiseStd*randn(p,1);
    y = rescale(t,0,1);
    X = zscore(X);
end

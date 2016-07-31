
function [startT, T, E, likelihood] = EM(X, m, n, itter)
    % X - 1 x k emission variables
    % m - ammount of possible states (y)ee
    % n - amount of possible emmissions (x)
    % initial estimation parameters
    % [startT, T, E] = genRandParams(m, n);
    % [T, E] = hmmtrain(X,T,E);
    % T(:, 1) = [];
    % startT = T(1, :);
    % T(1, :) = [];
    % E(1, :) = [];
    % return;
    
    k = length(X);
    epsilon = 10 ^ -4;
    [startT, T, E] = genRandParamsExt(m, n);
    % errors = [];
    likelihood = [];
    itters = 1:itter;
    for i=itters;

        % m x k
        [alpha, scale] = forwardAlg(X, startT, T, E);
        % m x k
        beta = backwardAlg(X, startT, T, E, scale);
        % m x k
        % gamma_t(i) = P(y_t = i|x_1:k)
        gamma = alpha .* beta ./ repmat(sum(alpha .* beta, 1), [m, 1]);
        
        xi = zeros(m);
        for t = 1 : k - 1
            newXi = (alpha(:, t) * (beta(:, t + 1) .* E(:, X(t + 1)))') .* T;
            xi = xi + newXi / sum(sum(newXi));
        end

        % update estimation parameters
        startT = gamma(:, 1);
        T =  xi;
        for i = 1 : n
            E(:, i) = E(:, i) + sum(gamma(:, X==i), 2);
        end
        T = bsxfun(@times, T, 1 ./ sum(T, 2));
        E = bsxfun(@times, E, 1 ./ sum(E, 2));
        likelihood(end + 1) = sum(log(scale));
        if length(likelihood)>1 & abs((likelihood(end) - likelihood(end -1)) / likelihood(end)) < epsilon
            % likelihood converged
            likelihood = likelihood(end);
            return
        end
    end
end



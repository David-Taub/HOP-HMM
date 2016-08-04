
function [startT, T, E, likelihood, gamma] = EM(X, m, n, itter);
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
    likelihood = -Inf;
    repeat = 3;
    for rep=1:repeat
        [repStartT, repT, repE] = genRandParamsExt(m, n);
        iterLike = [];
        for it=itters;

            % m x k
            [alpha, scale] = forwardAlg(X, repStartT, repT, repE);
            % m x k
            beta = backwardAlg(X, repStartT, repT, repE, scale);
            % m x k
            % gamma_t(i) = P(y_t = i|x_1:k)
            repGamma = alpha .* beta ./ repmat(sum(alpha .* beta, 1), [m, 1]);
            
            xi = zeros(m);
            for t = 1 : k - 1

                newXi = (alpha(:, t) * (beta(:, t + 1) .* repE(:, X(t + 1)))') .* repT;
                xi = xi + newXi / sum(sum(newXi));
            end

            % update estimation parameters
            repStartT = repGamma(:, 1);
            repT =  xi;
            for i = 1 : n
                repE(:, i) = repE(:, i) + sum(repGamma(:, X==i), 2);
            end
    	    repT = bsxfun(@times, repT, 1 ./ sum(repT, 2));
            repE = bsxfun(@times, repE, 1 ./ sum(repE, 2));
            iterLike(end + 1) = sum(log(scale));
            if length(iterLike)>1 & abs((iterLike(end) - iterLike(end -1)) / iterLike(end)) < epsilon
                % likelihood converged
                repLike = iterLike(end);
                fprintf('EM converged after %d iteratios: %f\n', it, repLike);
                break
            end
        end % end of itterations loop
        if likelihood < repLike
            likelihood = repLike;
            E = repE;
            T = repT;
            gamma = repGamma;
            startT = repStartT;
        end
    end % end of repeat loop
end



function [startT, T, E, likelihood, gamma] = EM(X, m, n, itter, tEpsilon)
    % X - 1 x k emission variables
    % m - ammount of possible states (y)ee
    % n - amount of possible emmissions (x)
    % initial estimation parameters
    
    k = length(X);
    epsilon = 10 ^ -4;
    bestLikelihood = -Inf;
    repeat = 3;
    for rep=1:repeat
        [startT, T, E] = genRandParamsExt(m, n);
        iterLike = [];
        for it=1:itter;

            % m x k
            [alpha, scale] = forwardAlg(X, startT, T, E);
            % m x k
            beta = backwardAlg(X, startT, T, E, scale);
            % m x k
            % gamma_t(i) = P(y_t = i|x_1:k)
            gamma = alpha .* beta ./ repmat(sum(alpha .* beta, 1), [m, 1]);
            
            xi = zeros(m);
            for t = 1 : k - 1
                % for each letter we get more information about the transition matrix update
                newXi = (alpha(:, t) * (beta(:, t + 1) .* E(:, X(t + 1)))') .* T;
                xi = xi + (newXi / sum(sum(newXi)));
            end

            % update estimation parameters
            startT = gamma(:, 1);
            T = T + xi;
            for i = 1 : n
                E(:, i) = E(:, i) + sum(gamma(:, X==i), 2);
            end
    	    T = bsxfun(@times, T, 1 ./ sum(T, 2));
            E = bsxfun(@times, E, 1 ./ sum(E, 2));
            for i = 1 : m
                for j = 1 : m
                    if i ~= j & T(i,j) > tEpsilon
                        T(i, i) = T(i, i) + (T(i, j) - tEpsilon);
                        T(i, j) = tEpsilon;
                    end
                end
            end
            iterLike(end + 1) = sum(log(scale));
            if length(iterLike)>1 & abs((iterLike(end) - iterLike(end -1)) / iterLike(end)) < epsilon
                % likelihood converged
                likelihood = iterLike(end);
                fprintf('EM converged after %d iteratios: %f\n', it, likelihood);
                break
            end
        end % end of itterations loop

        if bestLikelihood < likelihood
            bestLikelihood = likelihood;
            bestE = E;
            bestT = T;
            bestGamma = gamma;
            bestStartT = startT;
        end
    end % end of repeat loop
    
    % return parameters with best likelihood
    E = bestE;
    T = bestT;
    gamma = bestGamma;
    startT = bestStartT;
    likelihood = bestLikelihood;
end


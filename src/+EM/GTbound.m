
% G - m x k
% T - m x m
function [G, T, startT] = GTbound(params, G, T, doResample)
    if params.m == 1
        return;
    end
    G = exp(G);
    T = exp(T);
    originG = G;
    originT = T;
    if doResample
        G = makeDifferent(params, G);
    end

    fprintf('\n')
    for i=1:params.m
        fprintf('original: 1/1-T - %.2f 1/sum(G) - %.2f\n', 1./(1 - T(i, i)), 1./(sum(G(i, :), 2)))
    end
    % [G, T] = balanceGTweights(params, G, T);

    params.maxT < T
    params.minT > T
    sum(params.maxG < G, 2)
    sum(params.minG > G, 2)

    for i = 1:5
        T = max(T, params.minT);
        G = max(G, params.minG);
        T = min(T, params.maxT);
        G = min(G, params.maxG);
        s = sum(G, 2) + sum(T, 2);

        T = T ./ repmat(s, [1, params.m]);
        G = G ./ repmat(s, [1, params.k]);
    end
    % T = T + diag(1 - s);

    for i=1:params.m
        fprintf('balanceGTweights:  1/1-T - %.2f 1/sum(G) - %.2f\n', 1./(1 - T(i, i)), 1./(sum(G(i, :), 2)))
    end
    % T = limitTDiag(params, T);

    % for i=1:params.m
    %     fprintf('limitTDiag:  1/1-T - %.2f 1/sum(G) - %.2f\n', 1./(1 - T(i, i)), 1./(sum(G(i, :), 2)))
    % end
    T = limitCrossEnhancers(params, T, params.PCrossEnhancers);

    % for i=1:params.m
    %     fprintf('limitCrossEnhancers:  1/1-T - %.2f 1/sum(G) - %.2f\n', 1./(1 - T(i, i)), 1./(sum(G(i, :), 2)))
    % end

    startT = [ones(params.m-1, 1) * eps; 1 - (eps * (params.m - 1))];
    G = log(G);
    T = log(T);
    startT = log(startT);
end

function G = makeDifferent(params, G)

    [m, k] = size(G);
    for i = 1:m
        for j = i+1:m
            n1 = norm(G(i, :));
            n2 = norm(G(j, :));
            similarity = (G(i, :) ./ n1) * (G(j, :)./ n2)';
            if similarity > 0.95
                fprintf('state %d and %d are too similar (%.2f), resampling %d\n', i, j, similarity, j);
                G(j, :) = rand(1, k);
                G(j, :) = G(j, :) ./ (sum(G(j, :), 2) * params.PTotalBaseToSub);
            end
        end
    end
end
% T - m x m
function T = limitTDiag(params, T)
    for i = 1 : params.m - 1
        % limits enhancer length
        T = makeDiagonalDominant(params, T, i, params.m, params.PEnhancerToBackground);
        for j = 1 : params.m - 1
            % limits cross enhancers occurrences
            T = makeDiagonalDominant(params, T, i, j, params.PCrossEnhancers);
        end
    end
    for j = 1 : params.m - 1
        % limits background length
        T = makeDiagonalDominant(params, T, params.m, j, params.PBackgroundToEnhancer);
    end
end

% T - m x m
% transfers from T(i,j) to T(i,i)
function T = makeDiagonalDominant(params, T, i, j, threshold)
    if isExceedThreshold(params, threshold, T(i, j)) & i ~= j
        T(i, i) = T(i, i) + T(i, j) - threshold;
        T(i, j) = threshold;
    end
end


% T - m x m
% transfers from T(i,j) to T(i,i)
function T = limitCrossEnhancers(params, T, threshold)
    for i = 1:params.m - 1
        for j = 1:params.m - 1
            if j == i
                continue;
            end
            if T(i, j) > threshold
                T(i, params.m) = T(i, params.m) + T(i, j) - threshold;
                T(i, j) = threshold;
            end
        end
    end
end

% G - m x k
% T - m x m
function [G, T] = balanceGTweights(params, G, T)
    for i = 1:params.m - 1
        if isExceedThreshold(params, params.PTotalBaseToSub, sum(G(i, :), 2))
            G(i, :) = G(i, :) .* (params.PTotalBaseToSub / sum(G(i, :), 2));
            T(i, :) = T(i, :) .* ((1-params.PTotalBaseToSub) / (sum(T(i, :), 2) + eps));
        end
    end
    % make background with eps sub modes (motifs)
    G(params.m, :) = G(params.m, :) .* ((params.k*eps) / sum(G(params.m, :), 2));
    T(params.m, :) = T(params.m, :) .* ((1-(params.k*eps)) / (sum(T(params.m, :), 2) + eps));
end

function res = isExceedThreshold(params, thresh, val)
    res = abs(val - thresh) > thresh * params.maxPRatio;
end


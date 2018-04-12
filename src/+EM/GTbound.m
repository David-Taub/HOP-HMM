
% G - m x k
% T - m x m
function [G, T] = GTbound(params, G, T, doResample)
    if params.m == 1
        return;
    end
    G = exp(G);
    T = exp(T);
    originG = G;
    originT = T;
    % TODO: most of these bounding are probably not necessary, and could be removed
    if doResample
        G = makeDifferent(params, G);
    end
    [G, T] = balanceGTweights(params, G, T);
    fprintf('balanceGT, G Binding: %.3f T Binding: %.3f\n', mean(abs(originG(:) - G(:)), 1), mean(abs(originT(:) - T(:)), 1));
    T = limitTDiag(params, T);
    fprintf('limitTDiag, G Binding: %.3f T Binding: %.3f\n', mean(abs(originG(:) - G(:)), 1), mean(abs(originT(:) - T(:)), 1));
    G = log(G);
    T = log(T);
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
        T = transferWeight(params, T, i, params.m, params.PEnhancerToBackground);
        for j = 1 : params.m - 1
            T = transferWeight(params, T, i, j, params.PCrossEnhancers);
        end
    end
    for j = 1 : params.m - 1
        T = transferWeight(params, T, params.m, j, params.PBackgroundToEnhancer);
    end
end

% T - m x m
% transfers from T(i,j) to T(i,i)
function T = transferWeight(params, T, i, j, threshold)
    if isExceedThreshold(params, threshold, T(i, j)) & i ~= j
        T(i, i) = T(i, i) + T(i, j) - threshold;
        T(i, j) = threshold;
    end
end

% G - m x k
% T - m x m
function [G, T] = balanceGTweights(params, G, T)
    for i = 1:params.m-1
        if isExceedThreshold(params, params.PTotalBaseToSub, sum(G(i, :), 2))
            G(i, :) = G(i, :) .* (params.PTotalBaseToSub / sum(G(i, :), 2));
            T(i, :) = T(i, :) .* ((1-params.PTotalBaseToSub) / (sum(T(i, :), 2) + eps));
        end
    end
    G(params.m, :) = G(params.m, :) .* ((params.k*eps) / sum(G(params.m, :), 2));
    T(params.m, :) = T(params.m, :) .* ((1-(params.k*eps)) / (sum(T(params.m, :), 2) + eps));
end

function res = isExceedThreshold(params, thresh, val)
    res = abs(val - thresh) > thresh * params.maxPRatio;
end

% G - m x k
% T - m x m
function [G, T] = GTbound(params, G, T)
    if params.m == 1
        return;
    end
    G = exp(G);
    T = exp(T);
    [G, T] = balanceGTweights(params, G,T);
    T = limitTDiag(params, T);
    G = log(G);
    T = log(T);
end

% T - m x m
function T = limitTDiag(params, T)
    for i = 1 : params.m - 1
        T = transferWeight(T, i, params.m, params.maxPRatio*params.PEnhancerToBackground);
        for j = 1 : params.m - 1
            T = transferWeight(T, i, j, params.maxPRatio*params.PCrossEnhancers);
        end
    end
    for j = 1 : params.m - 1
        T = transferWeight(T, params.m, j, params.maxPRatio*params.PBackgroundToEnhancer);
    end
end

% T - m x m
function T = transferWeight(T, i, j, threshold)
    if T(i, j) > threshold & i ~= j
        T(i, i) = T(i, i) + T(i, j) - threshold;
        T(i, j) = threshold;
    end
end

% G - m x k
% T - m x m
function [G, T] = balanceGTweights(params, G, T)
    PTotalBaseToSub = misc.ratio2TransitionProb(mean(params.lengths, 2), params.enhancerMotifsRatio);
    for i = 1:params.m-1
        if sum(G(i, :), 2) > params.maxPRatio*params.PTotalBaseToSub
            G(i, :) = G(i, :) .* (params.PTotalBaseToSub / sum(G(i, :), 2));
            T(i, :) = T(i, :) .* ((1-params.PTotalBaseToSub) / (sum(T(i, :), 2) + eps));
        end
    end
    G(params.m, :) = G(params.m, :) .* ((params.k*eps) / sum(G(params.m, :), 2));
    T(params.m, :) = T(params.m, :) .* ((1-(params.k*eps)) / (sum(T(params.m, :), 2) + eps));
end
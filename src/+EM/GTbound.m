% help the algorithm converge to the correct result
% G - m x k
% T - m x m
function [G, T] = GTbound(params, G, T)
    if params.m == 1
        return;
    end
    assert(all(abs(sum(G, 2) + sum(T, 2) - 1) < 0.00001, 1))
    % G = exp(G);
    % T = exp(T);
    originG = G;
    originT = T;
    A = T; A(A > 0.5) = 1 - A(A > 0.5);round(1 ./ A)
    T(T > params.maxT) = params.maxT(T > params.maxT);
    T(T < params.minT) = params.minT(T < params.minT);
    G(G > params.maxG) = params.maxG(G > params.maxG);
    G(G < params.minG) = params.minG(G < params.minG);
    overG = sum(G, 2) > params.maxEnhMotif & [1:params.m]' <= params.enhancerAmount;
    G(overG, :) = params.maxEnhMotif * G(overG, :) ./ repmat(sum(G(overG, :), 2), [1, params.k]);
    underG = (sum(G, 2) < params.minEnhMotifTotal) & [1:params.m]' <= params.enhancerAmount;
    G(underG, :) = params.minEnhMotifTotal * G(underG, :) ./ repmat(sum(G(underG, :), 2), [1, params.k]);
    % bind sum of motif prob
    % T = T ./ repmat(sum(T, 2) + sum(G, 2), [1, params.m]);
    % G = G ./ repmat(sum(T, 2) + sum(G, 2), [1, params.k]);
    % T(eye(params.m) == 1) = T(eye(params.m) == 1) + 1 - sum(T, 2) - sum(G, 2);
    s = sum(params.minT, 2) + sum(params.minG, 2);
    Ta = T - params.minT;
    Ga = G - params.minG;
    % Ga = (params.maxEnhMotif .* (Ga + params.minG) ./ repmat(sum(Ga + params.minG, 2), [1, params.k])) - params.minG;
    sa = sum(Ta, 2) + sum(Ga, 2);
    Ta = (1 - s) .* Ta ./ repmat(sa, [1, params.m]);
    Ga = (1 - s) .* Ga ./ repmat(sa, [1, params.k]);
    T = params.minT + Ta;
    G = params.minG + Ga;


    A = T; A(A > 0.5) = 1 - A(A > 0.5);round(1 ./ A)

    assert(all(abs(sum(G, 2) + sum(T, 2) - 1) < 0.000001, 1));
    % assert(all(params.maxG(:) >= G(:) >= params.minG(:)));
    % assert(all(params.maxT(:) >= T(:) >= params.minT(:)));
    fprintf('%d/%d (G) and %d/%d (T) changed in bounding process\n', sum(sum(originG ~= G, 2), 1), params.m * params.k, sum(sum(originT ~= T, 2), 1), params.m ^ 2);
    G = log(G);
    T = log(T);

    assert(not(any(isnan(T(:)))));
    assert(not(any(isnan(G(:)))));
    assert(isreal(T));
    assert(isreal(G));
end


% T - m x m
% transfers from T(i,j) to T(i,i)
function T = limitCrossEnhancers(params, T, threshold)
    for i = 1:params.enhancerAmount
        for j = 1:params.enhancerAmount
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


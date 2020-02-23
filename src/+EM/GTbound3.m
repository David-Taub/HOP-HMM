% help the algorithm converge to the correct result
% G - m x k
% T - m x m
function [G, T, startT] = GTbound3(params, G, T, startT)
    EPS = 0.001;
    if params.m == 1
        return;
    end
    assert(all(abs(sum(G, 2) + sum(T, 2) - 1) < EPS, 1))
    % G = exp(G);
    % T = exp(T);
    originG = G;
    originT = T;
    originStartT = startT;

    MAX_ENH_LENGTH = 800;
    MIN_ENH_LENGTH = 20;
    % diag
    maxT = ones(params.m) * (1 / MIN_ENH_LENGTH);
    maxT(eye(params.m) == 1) = 1 - 1 / MAX_ENH_LENGTH;
    minT = eye(params.m) .* (1 - 1 / MIN_ENH_LENGTH);
    minT
    T
    maxT
    T(T > maxT) = maxT(T > maxT);
    T(T < minT) = minT(T < minT);

    s = sum(G, 2) + sum(T, 2);
    T = T ./ repmat(s, [1, params.m]);
    G = G ./ repmat(s, [1, params.k]);
    startT = startT ./ sum(startT, 1);
    assert(all(abs(sum(G, 2) + sum(T, 2) - 1) < EPS, 1));
    fprintf('%d/%d (G) and %d/%d (T) changed in bounding process\n', sum(sum(abs(originG - G) > EPS, 2), 1), ...
            params.m * params.k, sum(sum(abs(originT - T) > 0.001, 2), 1), ...
            params.m ^ 2);
    G = log(G);
    T = log(T);
    startT = log(startT);

    assert(not(any(isnan(T(:)))));
    assert(not(any(isnan(G(:)))));
    assert(not(any(isnan(startT(:)))));
    assert(isreal(T));
    assert(isreal(G));
    assert(isreal(startT));
end


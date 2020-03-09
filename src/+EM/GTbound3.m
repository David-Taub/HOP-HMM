% help the algorithm converge to the correct result
% G - m x k
% T - m x m

function theta = GTbound3(params, theta)

    for i = 1:10
        [G, T, startT] = GTbound4(params, exp(theta.G), exp(theta.T), exp(theta.startT));
    end
    theta.T = log(T);
    theta.G = log(G);
    theta.startT = log(startT);
end

function [G, T, startT] = GTbound4(params, G, T, startT)
    EPS = 0.001;
    if params.m == 1
        return;
    end
    originG = G;
    originT = T;
    originStartT = startT;

    startT(:) = 0;
    startT(end) = 1;
    T(T > params.maxT) = params.maxT(T > params.maxT);
    T(T < params.minT) = params.minT(T < params.minT);

    s = sum(G, 2) + sum(T, 2);
    T = T ./ repmat(s, [1, params.m]);
    G = G ./ repmat(s, [1, params.k]);
    startT = startT ./ sum(startT, 1);

    assert(all(abs(sum(G, 2) + sum(T, 2) - 1) < EPS, 1));
    fprintf('%d/%d (G) and %d/%d (T) changed in bounding process\n', sum(sum(abs(originG - G) > EPS, 2), 1), ...
            params.m * params.k, sum(sum(abs(originT - T) > 0.001, 2), 1), ...
            params.m ^ 2);

    assert(not(any(isnan(T(:)))));
    assert(not(any(isnan(G(:)))));
    assert(not(any(isnan(startT(:)))));
    assert(isreal(T));
    assert(isreal(G));
    assert(isreal(startT));
end


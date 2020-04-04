% help the algorithm converge to the correct result
% G - m x k
% T - m x m

function theta = GTbound3(params, theta)

    for i = 1:2
        [G, T, startT] = GTbound4(params, exp(theta.G), exp(theta.T), exp(theta.startT));
        % E = Ebound(params, exp(theta.E));
    end
    % theta.E = log(E);
    theta.T = log(T);
    theta.G = log(G);
    theta.startT = log(startT);
end

function E = Ebound(params, E)
    THRESHOLD = 1.3;
    for i = 1: params.enhancerAmount
        E(i, :) = max(E(1, :) ./ THRESHOLD, E(i, :));
        E(i, :) = min(E(1, :) .* THRESHOLD, E(i, :));
    end
end

function [G, T, startT] = GTbound4(params, G, T, startT)
    EPS = 0.00001;
    MIN_ENH_RATE = 130;
    if params.m == 1
        return;
    end
    originG = G;
    originT = T;
    originStartT = startT;
    MIN_G = [ones(params.enhancerAmount, 1) ./ MIN_ENH_RATE; ones(params.backgroundAmount, 1) .* EPS];
    MAX_G = [ones(params.enhancerAmount, 1); ones(params.backgroundAmount, 1) .* EPS];
    % Restrict startT
    startT(:) = EPS;
    startT(end) = 1 - EPS;

    % Restrict T
    T(T > params.maxT) = params.maxT(T > params.maxT);
    T(T < params.minT) = params.minT(T < params.minT);

    % Increase G
    TFFactor = ones(params.m, 1);
    TFFactor = max(TFFactor, MIN_G ./ sum(G, 2));
    TFFactor = min(TFFactor, MAX_G ./ sum(G, 2));
    G = G .* repmat(TFFactor, [1, params.k]);

    % L1 normalization
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


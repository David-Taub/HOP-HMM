function theta = genTheta(params, startTUniform, totalRandom)
    % normalized random probabilities
    if startTUniform || params.backgroundAmount == 0
        theta.startT = log(ones(params.m, 1) ./ params.m);
    else
        theta.startT = log([ones(params.enhancerAmount, 1) * params.EPS; ...
         ones(params.backgroundAmount, 1) * ...
         (1 - params.EPS * (params.enhancerAmount)) / params.backgroundAmount]);
    end
    if totalRandom
        M = rand(params.m, params.m + params.k);
        M = M ./ repmat(sum(M, 2), [1, params.m + params.k]);
        theta.T = M(:, 1:params.m);
        theta.G = M(:, params.m + 1:end);
        theta.T = log(theta.T);
        theta.G = log(theta.G);
    else
        theta.T = (params.maxT - params.minT) .* rand(params.m, params.m) + params.minT;
        theta.G = (params.maxG - params.minG) .* rand(params.m, params.k) + params.minG;

    end

    % m x n
    theta.E = rand([params.m, ones(1, params.order) .* params.n]);
    theta.E = log(bsxfun(@times, theta.E, 1 ./ sum(theta.E, params.order + 1)));

    theta = EM.GTbound3(params, theta);

    PRETRAINED_THETA_PATH = '../data/precomputation/pretrainedTheta.mat';
    if exist(PRETRAINED_THETA_PATH, 'file') == 2
        fprintf('Found pretrained theta file: %s\n', PRETRAINED_THETA_PATH);
        load(PRETRAINED_THETA_PATH);
        [foundM, foundK] = size(G);
        foundOrder = ndims(E) - 1;
        if params.m == foundM + 1 & params.k == foundK & foundOrder == params.order
            fprintf('Loading pretrained theta...\n');
            theta.E(1:foundM, :) = E(1:foundM, :);
            theta.G(1:foundM, 1:foundK) = G;
        end
    else
        fprintf('Using random theta initialization...\n');
    end

    assert(not(any(isnan(theta.T(:)))));
    assert(not(any(isnan(theta.E(:)))));
    assert(not(any(isnan(theta.G(:)))));
    assert(isreal(theta.T));
    assert(isreal(theta.E));
    assert(isreal(theta.G));
end
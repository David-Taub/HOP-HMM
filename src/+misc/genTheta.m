function [theta] = genTheta(params, isMixed)
    % normalized random probabilities
    if isMixed
        theta.startT = log(ones(params.m, 1) ./ params.m);
    else
        theta.startT = log([ones(params.m - 1, 1) * eps; 1 - eps * (params.m - 1)]);
    end
    theta.T = log((rand(params.m) .* (params.maxT - params.minT)) + params.minT);
    theta.G = log((rand(params.m, params.k) .* (params.maxG - params.minG)) + params.minG);

    % theta.startT = rand(params.m, 1);
    % theta.startT = ones(params.m, 1);
    % theta.startT = log(theta.startT / sum(theta.startT));



    % theta.T = eye(params.m) + params.PCrossEnhancers * rand(params.m);
    % theta.T = bsxfun(@times, theta.T, 1 ./ sum(theta.T, 2));

    % theta.G = rand(params.m, params.k);
    % theta.G = bsxfun(@times, theta.G, 1 ./ sum(theta.G, 2));
    % F = ones(params.m, 1) * 0.02;
    % theta.G = log(theta.G .* repmat(F, [1, params.k]));
    % theta.T = log(theta.T .* repmat(1-F, [1, params.m]));


    % m x n
    theta.E = rand([params.m, ones(1, params.order) .* params.n]);
    theta.E = log(bsxfun(@times, theta.E, 1 ./ sum(theta.E, params.order+1)));

    pretrainedThetaPath = '../data/precomputation/pretrainedTheta.mat';
    if exist(pretrainedThetaPath, 'file') == 2
        fprintf('Found pretrained theta file: %s\n', pretrainedThetaPath)
        load(pretrainedThetaPath);
        [foundM, foundK] = size(G);
        foundOrder = ndims(E) - 1;
        if params.m == foundM + 1 & params.k == foundK  & foundOrder == params.order
            fprintf('Loading pretrained theta...\n')
            theta.E(1:foundM, :) = E(1:foundM, :);
            theta.G(1:foundM, 1:foundK) = G;
        end
    else
        fprintf('Using random theta initialization...\n')
    end

    assert(not(any(isnan(theta.T(:)))))
    assert(not(any(isnan(theta.E(:)))))
    assert(not(any(isnan(theta.G(:)))))
end



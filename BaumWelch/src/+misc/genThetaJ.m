function [theta] = genThetaJ(params)
    % normalized random probabilities
    % theta.startT = rand(params.m, 1);
    theta.startT = ones(params.m, 1);
    theta.startT = log(theta.startT / sum(theta.startT));
    % m x m
    theta.T = eye(params.m) + params.tEpsilon * rand(params.m);
    theta.T = bsxfun(@times, theta.T, 1 ./ sum(theta.T, 2));

    theta.G = rand(params.m, params.k);
    theta.G = bsxfun(@times, theta.G, 1 ./ sum(theta.G, 2));
    F = ones(params.m, 1) * 0.02;
    theta.G = log(theta.G .* repmat(F, [1, params.k]));
    theta.T = log(theta.T .* repmat(1-F, [1, params.m]));


    % m x n
    theta.E = rand([params.m, ones(1, params.order) .* params.n]);
    theta.E = log(bsxfun(@times, theta.E, 1 ./ sum(theta.E, params.order+1)));
    % [params.PWMs, params.lengths] = misc.PWMs();
    % params.PWMs = log(params.PWMs);
    % params.lengths = lengths;
    % load(fullfile('data', 'dummyDNA.mat'));
end



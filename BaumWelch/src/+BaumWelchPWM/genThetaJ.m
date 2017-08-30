function [theta] = genThetaJ(params)
    % normalized random probabilities
    theta.startT = rand(params.m, 1);
    theta.startT = log(theta.startT / sum(theta.startT));
    % m x m
    theta.T = eye(params.m) + params.tEpsilon * rand(params.m);
    theta.T = log(bsxfun(@times, theta.T, 1 ./ sum(theta.T, 2)));

    theta.G = rand(params.m, params.k);
    theta.G = log(bsxfun(@times, theta.G, 1 ./ sum(theta.G, 2)));

    theta.F = log(ones(params.m, 1) * 0.03);
    % theta.F = log(ones(params.m, 1) * eps);

    % m x n
    theta.E = rand([params.m, ones(1, params.order) .* params.n]);
    theta.E = log(bsxfun(@times, theta.E, 1 ./ sum(theta.E, params.order+1)));
    [theta.PWMs, theta.lengths] = BaumWelchPWM.PWMs();
    theta.PWMs = log(theta.PWMs);
    % load(fullfile('data', 'dummyDNA.mat'));
end



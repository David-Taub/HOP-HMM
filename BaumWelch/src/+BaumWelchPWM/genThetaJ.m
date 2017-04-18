function [theta] = genThetaJ(params)
    % normalized random probabilities
    theta.startT = rand(params.m, 1);
    theta.startT = theta.startT / sum(theta.startT);
    % m x m
    theta.T = eye(params.m) + params.tEpsilon * rand(params.m);
    theta.T = bsxfun(@times, theta.T, 1 ./ sum(theta.T, 2));

    theta.M = rand(params.m, params.k);
    theta.M = bsxfun(@times, theta.M, 1 ./ sum(theta.M, 2));

    theta.F = rand(params.m, 1);

    % m x n
    theta.E = rand([params.m, ones(1, params.order) .* params.n]);
    theta.E = bsxfun(@times, theta.E, 1 ./ sum(theta.E, params.order+1));
    [theta.PWMs, theta.lengths] = BaumWelchPWM.PWMs();
end


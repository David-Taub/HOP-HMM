function [theta] = genThetaUni(params)
    % normalized random probabilities
    theta.startT = ones(params.m, 1);
    theta.startT = log(theta.startT / sum(theta.startT));
    % m x m
    theta.T = eye(params.m) + params.tEpsilon * ones(params.m);
    theta.T = bsxfun(@times, theta.T, 1 ./ sum(theta.T, 2));

    theta.G = ones(params.m, params.k);
    theta.G = bsxfun(@times, theta.G, 1 ./ sum(theta.G, 2));
    F = ones(params.m, 1) * 0.03;

    theta.G = log(theta.G .* repmat(F, [1, params.k]));
    theta.T = log(theta.T .* repmat(1-F, [1, params.m]));


    % m x n
    theta.E = ones([params.m, ones(1, params.order) .* params.n]);
    theta.E = log(bsxfun(@times, theta.E, 1 ./ sum(theta.E, params.order+1)));
end



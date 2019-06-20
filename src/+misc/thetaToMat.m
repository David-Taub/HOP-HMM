function mat = thetaToMat(params, theta, withT)
    mat = zeros(params.m, 1 + params.m + (params.n ^ params.order) + params.k);
    for i = 1:params.m
        mat(i, :) = [theta.T(i, :), theta.E(i, :), theta.G(i, :), theta.startT(i)];
    end
    if ~withT
        mat = mat(:, params.m + 1:end);
    end
end


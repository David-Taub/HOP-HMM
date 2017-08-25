

% alpha - N x m x L
% beta - N x m x L
% pX - N x 1
% xi - N x m x m x L
function xi = makeXi(theta, params, alpha, beta, X, pX)
    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % Eps - N x m x L
    Eps = BaumWelchPWM.EM.getEp3d(theta, params, X, 1:params.L, kronMN, matSize);
    compF = log(1-exp(theta.F));
    % xi - N x m x m x L
    xi = repmat(permute(alpha, [1, 2, 4, 3]), [1, 1, params.m, 1]);

    xi  = xi + repmat(compF', [params.N, 1, params.m, params.L]);
    xi  = xi + repmat(permute(theta.T, [3, 1, 2]), [params.N, 1, 1, params.L]);
    xi  = xi + repmat(permute(Eps, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    betaMoved = cat(3, beta(:, :, 2:end), -inf(params.N, params.m, 1));
    xi  = xi + repmat(permute(betaMoved, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    xi  = xi - repmat(pX, [1, params.m, params.m, params.L]);

end

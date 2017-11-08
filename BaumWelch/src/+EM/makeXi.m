

% alpha - N x m x L
% beta - N x m x L
% pX - N x 1
% xi - N x m x m x L
function xi = makeXi(theta, params, alpha, beta, X, pX)
    [N, L] = size(X);
    % Eps - N x m x L
    Eps = cat(3, EM.getEp3d(theta, params, X, 2:L), -inf(N, params.m, 1));
    % xi - N x m x m x L
    xi = repmat(permute(alpha, [1, 2, 4, 3]), [1, 1, params.m, 1]);

    xi  = xi + repmat(permute(theta.T, [3, 1, 2]), [N, 1, 1, L]);
    xi  = xi + repmat(permute(Eps, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    betaMoved = cat(3, beta(:, :, 2:end), -inf(N, params.m, 1));
    xi  = xi + repmat(permute(betaMoved, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    xi  = xi - repmat(pX, [1, params.m, params.m, L]);

end

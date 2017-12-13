% alpha - N x m x L
% beta - N x m x L
% X - N x L
% pcPWMp - N x k x L
% psi - N x m x k x L
function psi = makePsi(alpha, beta, X, params, theta, pcPWMp, pX)
    [N, L] = size(X);
    Eps = EM.getEp3d(theta, params, X, 1:L);
    Eps = cat(3, Eps, -inf(N, params.m, params.J + 1));
    pcPWMp = cat(3, pcPWMp, -inf(N, params.k, 1));
    beta = cat(3, beta, -inf(N, params.m, params.J + 1));
    % N x m x L
    % psi = zeros(N, params.m, params.k, L);
    psi = repmat(permute(theta.G, [3, 1, 2]), [N, 1, 1, L]);
    for t = 1:L
        % N x m x k x L - J - 1
        psi(:, :, :, t) = psi(:, :, :, t) + repmat(alpha(:, :, t), [1, 1, params.k]);
        % psi(:, :, :, t) = psi(:, :, :, t) + repmat(permute(theta.G, [3, 1, 2]), [N, 1, 1]);
        for l = 1:params.k
            % N x m x k
            psi(:, :, l, t) = psi(:, :, l, t) + repmat(pcPWMp(:, l, t+1), [1, params.m]);
            psi(:, :, l, t) = psi(:, :, l, t) + Eps(:, :, t+params.lengths(l)+1);
            psi(:, :, l, t) = psi(:, :, l, t) + beta(:, :, t+params.lengths(l)+1);
        end
    end
    psi  = psi - repmat(pX, [1, params.m, params.k, L]);
end

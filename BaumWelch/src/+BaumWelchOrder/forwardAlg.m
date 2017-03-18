
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
function [alpha, scale] = forwardAlg(Xs, startT, T, E)
    % alpha(N, i, j) P(y_s_j=i| x_s_1, ...x_s_j, startT, T, E)
    % scale(N, i) = P(x_s_i| startT, T, E)
    % m x L
    [N, L] = size(Xs);
    m = size(T, 1);
    kronMN = kron(1:m, ones(1, N));
    alpha = zeros(N, m, L);
    scale = zeros(N, L);
    order = matDim(E) - 1;
    startE = sumDim(E, [2 : order]);
    alpha(:, :, 1) = (repmat(startT, [1, N]) .* startE(:, Xs(:, 1))).';
    scale(:, 1) = sum(alpha(:, :, 1), 2);
    matSize = [m , 4 * ones(1, order)];
    for t = 2 : L
        % N x m
        Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % N x m * m x m .* N x m = N x m
        newAlphas = (alpha(:, :, t-1) * T) .* Ep;
        % N x 1
        scale(:, t) = sum(newAlphas, 2);
        alpha(:, :, t) = bsxfun(@times, newAlphas, 1 ./ scale(:, t));
    end
end

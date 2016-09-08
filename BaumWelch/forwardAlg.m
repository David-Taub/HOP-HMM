
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - S x k emission variables
function [alpha, scale] = forwardAlg(Xs, startT, T, E)
    % alpha(s, i, j) P(y_s_j=i| x_s_1, ...x_s_j, startT, T, E)
    % scale(s, i) = P(x_s_i| startT, T, E)
    % m x k
    [S, k] = size(Xs);
    m = length(T);
    alpha = zeros(S, m, k);
    scale = zeros(S, k);
    order = matDim(E) - 1;
    startE = sumDim(E, [2 : order]);
    alpha(:, :, 1) = (repmat(startT, [1, S]) .* startE(:, Xs(:, 1))).';
    scale(:, 1) = sum(alpha(:, :, 1), 2);
    matSize = [m , 4 * ones(1, order)];
    for t = 2 : k
        % newAlphas = E(:, X(t)) .* sum(repmat(alpha(:, t-1), [1, m]) .* T, 1).';
        Et = E;
        if t < order
            Et = sumDim(E, [2 : 1 + order - t]);
            subscripts = [kron(1:m, ones(1, S)); repmat(Xs(:, 1 : t).', [1, m])];
            indices = matSub2ind(matSize(1 : t + 1), subscripts);
        else
            subscripts = [kron(1:m, ones(1, S)); repmat(Xs(:, t-order+1 : t).', [1, m])];
            indices = matSub2ind(matSize, subscripts);
        end

        % newAlphas = (T.' * alpha(:, t-1)) .* Et(indices).';
        % scale(t) = sum(newAlphas);
        % alpha(:, t) = newAlphas / scale(t);
        % S x m
        Ep = reshape(Et(indices).', [S, m]);
        % S x m * m x m .* S x m = S x m
        newAlphas = (alpha(:, :, t-1) * T) .* Ep;
        % S x 1
        scale(:, t) = sum(newAlphas, 2);
        alpha(:, :, t) = bsxfun(@times, newAlphas, 1 ./ scale(:, t));
    end
end

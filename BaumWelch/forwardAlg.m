
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - 1 x k emission variables
function [alpha, scale] = forwardAlg(X, startT, T, E)
    % alpha(i,j) P(y_j=i| x_1, ...x_j, startT, T, E)
    % scale(i) = P(x_i| startT, T, E)
    % meaning: alpha_j(i) = alpha(j, i)
    % m x k
    k = length(X);
    m = length(T);
    alpha = zeros(m, k);
    scale = zeros(1, k);
    order = matDim(E) - 1;
    startE = sumDim(E, [2 : order]);
    alpha(:, 1) = startT .* startE(:, X(1));
    scale(1) = sum(alpha(:, 1));
    matSize = [m , 4 * ones(1, order)];
    for t = 2 : k
        % newAlphas = E(:, X(t)) .* sum(repmat(alpha(:, t-1), [1, m]) .* T, 1).';
        Et = E;
        if t < order
            Et = sumDim(E, [2 : 1 + order - t]);
            subscripts = [1:m ; repmat(X(1 : t).', [1, m])];
            indices = matSub2ind(matSize(1 : t + 1), subscripts);
        else
            subscripts = [1:m ; repmat(X(t-order+1 : t).', [1, m])];
            indices = matSub2ind(matSize, subscripts);
        end
        newAlphas = (T.' * alpha(:, t-1)) .* Et(indices).';
        scale(t) = sum(newAlphas);
        alpha(:, t) = newAlphas / scale(t);
    end
    % checked
end

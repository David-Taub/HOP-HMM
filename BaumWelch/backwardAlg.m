
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - k x 1 emission variables
function beta = backwardAlg(X, startT, T, E, scale, eta)

    % beta(i,t) P( x_t+1, ...x_k| y_t=i, startT, T, E)

    k = length(X);
    m = length(T);
    % m x k
    order = matDim(E) - 1;
    matSize = [m , 4 * ones(1, order)];
    beta = ones(m, k) / m;

    for t = k - 1 : -1 : 1
        % newBetas = sum(repmat(E(:, X(t+1)) .* beta(:, 1), [1, m]).' .* T, 2);
        % beta= cat(2, newBetas, beta);
        Et = E;
        if t < order
            Et = sumDim(E, [2 : 1 + order - t]);
            subscripts = [1:m ; repmat(X(1 : t).', [1, m])];
            indices = matSub2ind(matSize(1 : t + 1), subscripts);
        else
            subscripts = [1:m ; repmat(X(t-order+1 : t).', [1, m])];
            indices = matSub2ind(matSize, subscripts);
        end
        beta(:, t) = (T * Et(indices).' .* (beta(:, t+1) .^ eta)) / scale(t);
        % beta(:, t) = (T * (E(:, X(t + 1)) .* beta(:, t+1))) / scale(t);
    end
    % checked
end



% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - S x k emission variables
% beta(s, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)

function beta = backwardAlg(Xs, startT, T, E, scale)
    [S, k] = size(Xs);
    m = length(T);
    % m x k
    order = matDim(E) - 1;
    matSize = [m , 4 * ones(1, order)];
    beta = ones(S, m, k);

    for t = k - 1 : -1 : 1
        % newBetas = sum(repmat(E(:, X(t+1)) .* beta(:, 1), [1, m]).' .* T, 2);
        % beta= cat(2, newBetas, beta);
        Et = E;
        if t < order
            Et = sumDim(E, [2 : 1 + order - t]);

            subscripts = [kron(1:m, ones(1, S)); repmat(Xs(:, 1 : t).', [1, m])];
            indices = matSub2ind(matSize(1 : t + 1), subscripts);
        else
            subscripts = [kron(1:m, ones(1, S)); repmat(Xs(:, t-order+1 : t).', [1, m])];
            indices = matSub2ind(matSize, subscripts);
        end
        % S x m
        Ep = reshape(Et(indices).', [S, m]);
        % S x m 
        beta(:, :, t) = bsxfun(@times, (Ep .* beta(:, :, t+1)) * T.', 1 ./ scale(:, t));

    end
    % checked
end


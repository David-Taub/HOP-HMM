
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)

function beta = backwardAlg(Xs, startT, T, E, scale)
    [N, L] = size(Xs);
    m = length(T);
    kronMN = kron(1:m, ones(1, N));
    n = 4;
    % m x L
    order = matDim(E) - 1;
    matSize = [m , n * ones(1, order)];
    beta = ones(N, m, L);

    for t = L : -1 : 2
        Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % N x m 
        beta(:, :, t-1) = bsxfun(@times, (Ep .* beta(:, :, t)) * T.', 1 ./ scale(:, t-1));

    end
    % checked
end


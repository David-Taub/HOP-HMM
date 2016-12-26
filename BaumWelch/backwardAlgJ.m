
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(Xs, startT, T, E, scale, PWMsRep, lengths)
    [N, L] = size(Xs);
    m = length(T);
    [~, J, n, ~] = size(PWMsRep);
    % m x L
    order = matDim(E) - 1;
    matSize = [m , n * ones(1, order)];
    beta = ones(N, m, L);

    Xs1H = cat(2, mat23Dmat(Xs, n), zeros(N, J, n));
    for t = L : -1 : 2
        
        % N x m
        Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % N x m 
        newBeta = (Ep .* beta(:, :, t+1)) * T.';

        % beta is built by the base modes emission and the PWM submodes emission
        % N x m x k
        betaSlice = beta(:, :, t - 1 + lengths);
        newBeta = newBeta + PWMstep(betaSlice, PWMsRep, Xs1H, Y, J, t);
        beta(:, :, t-1) = bsxfun(@times, newBeta, 1 ./ scale(:, t-1));
        
    end
end


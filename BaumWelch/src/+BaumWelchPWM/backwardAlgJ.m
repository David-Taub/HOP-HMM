
% T - m x m transfer probability matrix between mode bases
% Y - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(Xs, startT, T, Y, F, E, scale, PWMsRep, lengths)
    [N, L] = size(Xs);
    [m, k] = length(Y);
    [~, J, n, ~] = size(PWMsRep);
    % m x L
    order = matDim(E) - 1;
    matSize = [m , n * ones(1, order)];
    beta = ones(N, m, L);
    Y = bsxfun(@times, Y, F);
    E = bsxfun(@times, E, 1-F);
    Xs1H = cat(2, mat23Dmat(Xs, n), zeros(N, J, n));

    for t = L : -1 : 2
        % N x m
        % note: this looks at part of the sequences before t, which might be problematic.
        % TODO: I should remove this note if everything goes well
        Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % N x m 
        newBeta = (Ep .* beta(:, :, t)) * T.';

        % beta is built by the base modes emission and the PWM submodes emission
        % N x m x k
        betaSlice = beta(:, :, t + lengths);
        newBeta = newBeta + PWMstep(betaSlice, PWMsRep, Xs1H, Y, J, t-1);
        beta(:, :, t-1) = bsxfun(@times, newBeta, 1 ./ scale(:, t-1));
        
    end
end


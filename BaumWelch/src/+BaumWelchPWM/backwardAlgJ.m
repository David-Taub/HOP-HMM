
% T - m x m transfer probability matrix between mode bases
% Y - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(Xs, T, Y, F, E, scale, lengths, pcPWMp, J)
    [N, L] = size(Xs);
    [m, k] = size(Y);
    % m x L
    order = matUtils.matDim(E) - 1;
    n = size(E, order + 1);

    kronMN = kron(1:m, ones(1, N));
    matSize = [m , n * ones(1, order)];
    beta = cat(3, ones(N, m, L), zeros(N, m, J));
    Y = bsxfun(@times, Y, F);
    E = bsxfun(@times, E, 1-F);
    % performance optimization
    Ys = shiftdim(Y, -1);
    Ys = repmat(Ys, [N, 1, 1]);

    for t = L : -1 : 2
        fprintf('Backward algorithm %.2f%%\r', 100 * (L - t + 1) / L);
        % N x m
        % note: this looks at part of the sequences before t, which might be problematic.
        % TODO: I should remove this note if everything goes well
        Ep = BaumWelchPWM.getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % N x m
        newBeta = (Ep .* beta(:, :, t)) * T.';

        % beta is built by the base modes emission and the PWM submodes emission
        % N x m x k
        betaSlice = beta(:, :, t + lengths - 1);

        % N x m
        newBeta = newBeta + BaumWelchPWM.PWMstep(betaSlice, Ys, repmat(t-1, [k, 1]), pcPWMp, J);
        % m x L
        beta(:, :, t-1) = bsxfun(@times, newBeta, 1 ./ scale(:, t-1));
    end
    beta = beta(:, :, 1:L);
    fprintf('\n')
end


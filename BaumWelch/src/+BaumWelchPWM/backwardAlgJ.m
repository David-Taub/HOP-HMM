
% T - m x m transfer probability matrix between mode bases
% theta.M - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% theta.E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(Xs, theta, params, scale, pcPWMp)

    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    beta = cat(3, ones(params.N, params.m, params.L), zeros(params.N, params.m, params.J));
    theta.M = bsxfun(@times, theta.M, theta.F);
    theta.E = bsxfun(@times, theta.E, 1-theta.F);
    % performance optimization
    Ms = shiftdim(theta.M, -1);
    Ms = repmat(Ms, [params.N, 1, 1]);

    for t = params.L : -1 : 2
        fprintf('Backward algorithm %.2f%%\r', 100 * (params.L - t + 2) / params.L);
        % N x m
        % note: this looks at part of the sequences before t, which might be problematic.
        % TODO: I should remove this note if everything goes well
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);
        % N x m
        newBeta = (Ep .* beta(:, :, t)) * theta.T.';

        % beta is built by the base modes emission and the PWM submodes emission
        % N x m x k
        betaSlice = beta(:, :, t + theta.lengths - 1);

        % N x m
        PWMstepP = BaumWelchPWM.PWMstep(betaSlice, Ms, repmat(t, [params.k, 1]), pcPWMp, params.J);
        newBeta = newBeta + PWMstepP * theta.T.';
        % m x L
        beta(:, :, t-1) = bsxfun(@times, newBeta, 1 ./ scale(:, t-1));

    end
    fprintf('\n')
end



% T - m x m transfer probability matrix between mode bases
% theta.G - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% theta.E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - N x L emission variables
% beta - N x m x L
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(X, theta, params, pcPWMp)

    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % zero appended to handle pwm steps in the end of the sequence (first iterations) which have probability 0
    beta = cat(3, zeros(params.N, params.m, params.L), -inf(params.N, params.m, params.J));
    % performance optimization
    Gs = repmat(shiftdim(theta.G, -1), [params.N, 1, 1]);
    Fs = repmat(theta.F.', [params.N, 1, params.k]);
    compF = repmat(log(1-exp(theta.F))', [params.N, 1]);
    expT = exp(theta.T.');
    for t = params.L : -1 : 2
        % fprintf('Backward algorithm %.2f%%\r', 100 * (params.L-t+2) / params.L);
        % note: this peeks at part of the sequences before t, which might be problematic
        % note 25.07.17: Tommy thinks it is fine - and I see no reason it will affect non-margins areas.
        % N x m
        Ep = BaumWelchPWM.EM.getEp(theta, params, X, t, kronMN, matSize);

        % No  1 2 PWM2 PWM2 2 1
        % Yes 1 2 PWM2 2 1 <== in non all to all mode, only these hidden states is legal
        % No  1 2 PWM2 1 1
        % No  1 1 PWM2 2 1
        % No  1 1 PWM2 1 1
        % No  1 1 PWM2 PWM2 1 1

        % N x m x k
        EpReturn = BaumWelchPWM.EM.getEp3d(theta, params, X, t+theta.lengths, kronMN, matSize);
        % N x m x k
        betaSlice = beta(:, :, t+theta.lengths);
        % N x m
        subStateStep = BaumWelchPWM.EM.PWMstep(betaSlice, Gs, repmat(t, [params.k, 1]), pcPWMp, EpReturn, Fs);
        % N x m
        baseStateStep = matUtils.logMatProd(Ep + beta(:, :, t), expT') + compF;
        beta(:,:,t-1) = matUtils.logAdd(baseStateStep, subStateStep);
    end
    beta = beta(:, :, 1:params.L);
    fprintf('\n');
end
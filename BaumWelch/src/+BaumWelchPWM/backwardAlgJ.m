
% T - m x m transfer probability matrix between mode bases
% theta.G - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% theta.E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta - N x m x L
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(Xs, theta, params, pcPWMp)

    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % zero appended to handle pwm steps in the end of the sequence (first iterations) which have probability 0
    beta = cat(3, zero(params.N, params.m, params.L), -inf(params.N, params.m, params.J));
    % performance optimization
    Gs = repmat(shiftdim(theta.G, -1), [params.N, 1, 1]);
    compF = log(1-exp(theta.F));
    expT = exp(theta.T.');
    for t = params.L : -1 : 2
        fprintf('Backward algorithm %.2f%%\r', 100 * (params.L-t+2) / params.L);
        % N x m
        % note: this peeks at part of the sequences before t, which might be problematic
        % note 25.07.17: Tommy thinks it is fine - and I see no reason it will affect non-margins areas.
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);

        % No  1 2 PWM2 PWM2 2 1
        % Yes 1 2 PWM2 2 1 <== in non all to all mode, only these hidden states is legal
        % No  1 2 PWM2 1 1
        % No  1 1 PWM2 2 1
        % No  1 1 PWM2 1 1
        % No  1 1 PWM2 PWM2 1 1

        % N x m x k
        EpReturn = BaumWelchPWM.getEp3d(theta, params, Xs, t+params.lengths, kronMN, matSize);
        % N x m x k
        betaSlice = beta(:, :, t+params.lengths);
        % N x m
        subStateStep = BaumWelchPWM.PWMstep(betaSlice, Gs, repmat(t, [params.k, 1]), pcPWMp, EpReturn, theta.F);
        baseStateStep = log(exp(Ep + beta(:, :, t)) * expT') + compF;
        beta(:,:,t-1) = matUtils.logAdd(baseStateStep, subStateStep);
    end
    beta = beta(:, :, 1:params.L);
    fprintf('\n');
end
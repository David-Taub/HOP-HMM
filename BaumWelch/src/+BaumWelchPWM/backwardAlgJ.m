
% T - m x m transfer probability matrix between mode bases
% theta.M - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% theta.E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta - N x m x L
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(Xs, theta, params, scale, pcPWMp)

    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % zero appended to handle pwm steps in the end of the sequence (first iterations) which have probability 0
    beta = cat(3, ones(params.N, params.m, params.L), zeros(params.N, params.m, params.J));
    if params.all2allMode
        theta.M = bsxfun(@times, theta.M, theta.F);
        theta.E = bsxfun(@times, theta.E, 1-theta.F);
    end
    % performance optimization
    Ms = shiftdim(theta.M, -1);
    Ms = repmat(Ms, [params.N, 1, 1]);
    for t = params.L : -1 : 2
        fprintf('Backward algorithm %.2f%%\r', 100 * (params.L-t+2) / params.L);
        % we are building the t-1 line, using t
        % N x m
        % note: this peeks at part of the sequences before t, which might be problematic
        % note 25.07.17: Tommy thinks it is fine - and I see no reason it will affect non-margins areas.
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);
        if params.all2allMode
            % N x m x k
            betaSlice = beta(:,:,t+params.lengths-1);
            % N x m x k
            PWMstepP = BaumWelchPWM.PWMstep(betaSlice, Ms, repmat(t, [params.k, 1]), pcPWMp);
            % N x m x k
            % Yes 1 2 PWM2 PWM2 2 1
            % Yes 1 2 PWM2 2 1
            % Yes 1 2 PWM2 1 1
            % Yes 1 1 PWM2 2 1
            % Yes 1 1 PWM2 1 1
            % Yes 1 1 PWM2 PWM2 1 1
            newBeta = (Ep .* beta(:,:,t) + sum(PWMstepP, 3)) * theta.T.';
        else
            % No  1 2 PWM2 PWM2 2 1
            % Yes 1 2 PWM2 2 1 <== in non all to all mode, only these hidden states is legal
            % No  1 2 PWM2 1 1
            % No  1 1 PWM2 2 1
            % No  1 1 PWM2 1 1
            % No  1 1 PWM2 PWM2 1 1

            % N x m x k
            EpReturn = BaumWelchPWM.getEp3d(theta, params, Xs, t+lengths, kronMN, matSize);
            % N x m x k
            betaSlice = beta(:,:,t+params.lengths);
            % N x m x k
            PWMstepP = BaumWelchPWM.PWMstep(betaSlice, Ms, repmat(t, [params.k, 1]), pcPWMp);
            Tf = theta.T;
            Tf(eye(m)==1) = Tf(eye(m)==1) .* (1-theta.F);
            % N x m = (N x m .* N x m) * (m x m) + (N x m .* N x m)
            newBeta = (Ep .* beta(:,:,t)) * Tf.' + sum(PWMstepP .* EpReturn, 3) .* repmat((diag(theta.T) .* theta.F)', [N, 1]);
        end
        beta(:,:,t-1) = bsxfun(@times, newBeta, 1 ./ scale(:,t-1));
    end
    beta = beta(:, :, 1:params.L);
    fprintf('\n');
end

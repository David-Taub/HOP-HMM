
% T - m x m transfer probability matrix between mode bases
% theta.M - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% theta.E - m x n emission matrix E_ij means x_t = j | y_t = i
% Xs - N x L emission variables
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(Xs, theta, params, scale, pcPWMp)

    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    betaBase = ones(params.N, params.m, params.L);
    betaSub = cat(3, ones(params.N, params.m, params.L), zeros(params.N, params.m, params.J));
    if params.all2allMode
        theta.M = bsxfun(@times, theta.M, theta.F);
        theta.E = bsxfun(@times, theta.E, 1-theta.F);
    end
    % performance optimization
    Ms = shiftdim(theta.M, -1);
    Ms = repmat(Ms, [params.N, 1, 1]);

    for t = params.L : -1 : 2
        fprintf('Backward algorithm %.2f%%\r', 100 * (params.L-t+2) / params.L);
        % N x m
        % note: this peeks at part of the sequences before t, which might be problematic.
        % TODO: I should remove this note if everything goes well
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);
        % N x m x k
        betaSlice = betaSub(:,:,sub2ind([params.L, params.k], t+lengths-1, [1:params.k]));
        % N x m x k
        PWMstepP = BaumWelchPWM.PWMstep(betaSlice, Ms, repmat(t, [params.k, 1]), pcPWMp, params.J);
        if params.all2allMode
            % N x m x k
            newBetaBase = (Ep .* betaBase(:,:,t) + sum(PWMstepP, 3)) * theta.T.';
        else
            Tf = theta.T;
            Tf(eye(m)==1) = Tf(eye(m)==1) .* (1-theta.F);
            newBetaBase = (Ep .* betaBase(:,:,t)) * Tf.' + sum(PWMstepP, 3) .* repmat((diag(theta.T) .* theta.F)', [N, 1]);
        end

        % beta is built by the base modes emission and the PWM submodes emission
        % N x m x k
        if params.all2allMode
            betaSlice = betaBase(:,:,t+theta.lengths-1) + sum(betaSub(:,:,t+theta.lengths-1,:), 4)
            beta(:,:,t+theta.lengths-1);
            % #########################################################
            % TODO: problem here, are all submodes equal?
            newBetaBase = (Ep .* betaBase(:,:,t) + sum(PWMstepP, 3)) * theta.T.';
        else
            betaSlice = betaBase(:,:,t+theta.lengths-1) + sum(betaSub(:,:,t+theta.lengths-1,:), 4)
            
        end
        1
        % N x m
        PWMstepP = BaumWelchPWM.PWMstep(betaSlice, Ms, repmat(t, [params.k, 1]), pcPWMp, params.J);
        newBeta = newBeta + PWMstepP * theta.T.';
        % m x L
        beta(:,:,t-1) = bsxfun(@times, newBeta, 1 ./ scale(:,t-1));

    end
    fprintf('\n')
end


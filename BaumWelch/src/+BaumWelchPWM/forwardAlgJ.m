% This model assumes m modes bases that emits like high order HHM,
% where each base can transfer into k submodes that emits with PWM from JASPAR project
% m - sum of enhancer and background modes (not a parameter)
% T - m x m transfer probability matrix between mode bases
% theta.G - m x k transfer probability matrix between mode bases and their PWM modes
% theta.F - m x 1 probability to get into a PWM mode per state.
% n - number of alphabet (4, i.e. ACGT)
% startT - m x 1 probabilities of first states
% theta.E - m x n x n x ... x n (order times) alphabet emission probability matrix
% PWMsRep - N x J x n x k: N replication of emission matrix of m Jaspar PWMs with length J
%        true length of i'th PWM< J and given in lengths(i) if a PWM is
%        shorter than j, it is aligned to the end of the 3rd dimension.
% lengths - m x 1 length of each motif in the PWM matrix. J = max(lengths)
% Xs - N x L emission variables
% alpha - N x m x L
function alpha = forwardAlgJ(Xs, theta, params, pcPWMp)
    matSize = [params.m , params.n * ones(1, params.order)];
    % m x L
    kronMN = kron(1:params.m, ones(1, params.N));
    compF = repmat(log(1-exp(theta.F))', [params.N, 1]);
    Gs = repmat(shiftdim(theta.G, -1), [params.N, 1, 1]);
    Fs = repmat(theta.F.', [params.N, 1, params.k]);
    expT = exp(theta.T);
    % the k+1 index is for base modes, 1 to k are for sub modes
    alpha = -inf(params.N, params.m, params.L + params.J);
    % N x m
    Ep = BaumWelchPWM.getEp(theta, params, Xs, 1, kronMN, matSize);
    alpha(:, :, params.J+1) = (repmat(theta.startT', [params.N, 1]) + Ep);

    % if length is 3, J = 4
    % then B is the base mode, S is the submode (PWM mode)
    % 654321t     - indices
    % BBBBSSS???? - hidden
    % XXXX123???? - emission
    for t = 2:params.L
        fprintf('Forward algorithm %.2f%%\r', 100*t/params.L);
        % N x m
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);
        % N x m
        baseStateStep = log(exp(alpha(:, :, params.J+t-1) + compF) * expT) + Ep;
        % N x m x k
        alphaSlice = alpha(:, :, params.J+t-theta.lengths-1);
        % N x m x k
        subStateStep = BaumWelchPWM.PWMstep(alphaSlice, Gs, t-theta.lengths', pcPWMp, repmat(Ep, [1, 1, params.k]), Fs);
        % N x m x k -> N x m
        subStateStep = matUtils.logMatSum(subStateStep, 3);
        alpha(:, :, params.J + t) = matUtils.logAdd(baseStateStep, subStateStep);
    end
    fprintf('\n');
    alpha = alpha(:, :, params.J+1:end);
end


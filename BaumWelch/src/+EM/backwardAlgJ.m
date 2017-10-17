
% T - m x m transfer probability matrix between mode bases
% theta.G - m x k transfer probability matrix between mode bases and their PWM modes
% startT - m x 1 probabilities of first states
% theta.E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - N x L emission variables
% pcPWMp - N x k x L
% beta - N x m x L
% beta(N, i, t) P( x_s_t+1, ...x_s_k| y_s_t=i, startT, T, E)
function beta = backwardAlgJ(X, theta, params, pcPWMp)
    [N, L] = size(X);
    kronMN = kron(1:params.m, ones(1, N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % zero appended to handle pwm steps in the end of the sequence (first iterations) which have probability 0
    beta = cat(3, zeros(N, params.m, L), -inf(N, params.m, params.J));
    % performance optimization
    Gs = repmat(shiftdim(theta.G, -1), [N, 1, 1]);
    expT = exp(theta.T.');
    for t = L : -1 : 2
        % fprintf('Backward algorithm %.2f%%\r', 100 * (L-t+2) / L);
        % note: this peeks at part of the sequences before t, which might be problematic
        % note 25.07.17: Tommy thinks it is fine - and I see no reason it will affect non-margins areas.
        % N x m
        Ep = EM.getEp(theta, params, X, t, kronMN, matSize);


        % N x m x k
        EpReturn = EM.getEp3d(theta, params, X, t+theta.lengths, kronMN, matSize);
        % N x m x k
        betaSlice = beta(:, :, t+theta.lengths);
        % N x m
        subStateStep = EM.PWMstep(betaSlice, Gs, repmat(t, [params.k, 1]), pcPWMp, EpReturn);
        % N x m
        baseStateStep = matUtils.logMatProd(Ep + beta(:, :, t), theta.T');
        beta(:,:,t-1) = matUtils.logAdd(baseStateStep, subStateStep);
    end
    beta = beta(:, :, 1:L);
    % fprintf('\n');
end
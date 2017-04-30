% This model assumes m modes bases that emits like high order HHM,
% where each base can transfer into k submodes that emits with PWM from JASPAR project
% m - sum of enhancer and background modes (not a parameter)
% T - m x m transfer probability matrix between mode bases
% theta.M - m x k transfer probability matrix between mode bases and their PWM modes
% theta.F - m x 1 probability to get into a PWM mode per state.
% n - number of alphabet (4, i.e. ACGT)
% startT - m x 1 probabilities of first states
% theta.E - m x n x n x ... x n (order times) alphabet emission probability matrix
% PWMsRep - N x J x n x k: N replication of emission matrix of m Jaspar PWMs with length J
%        true length of i'th PWM< J and given in lengths(i) if a PWM is
%        shorter than j, it is aligned to the end of the 3rd dimension.
% lengths - m x 1 length of each motif in the PWM matrix. J = max(lengths)
% Xs - N x L emission variables

function [alpha, alphaT, scale] = forwardAlgJ(Xs, theta, params, pcPWMp)
    % alpha(N, i, j) P(y_s_j=i| x_s_1, ...x_s_j, startT, T, PWMs)
    % scale(N, i) = P(x_s_i| startT, T, PWMs)
    % m x L

    kronMN = kron(1:params.m, ones(1, params.N));
    alpha = zeros(params.N, params.m, params.L + params.J);
    alphaT = zeros(params.N, params.m, params.L + params.J);
    scale = zeros(params.N, params.L);
    startE = matUtils.sumDim(theta.E, 2 : params.order);

    alpha(:, :, params.J+1) = (repmat(theta.startT, [1, params.N]) .* startE(:, Xs(:, 1))).';
    alphaT(:, :, params.J+1) = alpha(:, :, params.J+1) * theta.T;
    scale(:, 1) = sum(alpha(:, :, params.J+1), 2);
    theta.M = bsxfun(@times, theta.M, theta.F);
    theta.E = bsxfun(@times, theta.E, 1-theta.F);

    % performance optimization
    Ms = shiftdim(theta.M, -1);
    Ms = repmat(Ms, [params.N, 1, 1]);
    matSize = [params.m , params.n * ones(1, params.order)];
    % if length is 3, J = 4
    % then B is the base mode, S is the submode (PWM mode)
    % 654321t     - indices
    % BBBBSSS???? - hidden
    % XXXX123???? - emission
    for t = 2:4%params.L
        assert(not(any(isnan(alpha(:)))))
        fprintf('Forward algorithm %.2f%%\r', 100 * t / params.L);
        % N x m
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);
        % N x m
        newAlphas = alphaT(:, :, t + params.J  - 1) .* Ep;
        % N x m x k
        alphaSlice = alphaT(:, :, t + params.J - theta.lengths);
        PWMp = BaumWelchPWM.PWMstep(alphaSlice, Ms, t - theta.lengths' + 1, pcPWMp, params.J);
        newAlphas = newAlphas + PWMp;

        % N x 1
        scale(:, t) = sum(newAlphas, 2);
        alpha(:, :, params.J + t) = bsxfun(@times, newAlphas, 1 ./ scale(:, t));
        alphaT(:, :, params.J + t) = alpha(:, :, params.J + t) * theta.T;
    end
    fprintf('\n');

    alpha = alpha(:, :, params.J+1:end);
    alpha = cat(3, alpha, zeros(params.N, params.m, params.J));
end


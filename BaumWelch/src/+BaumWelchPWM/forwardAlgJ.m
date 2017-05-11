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

function [alphaBase, alphaSub, scale] = forwardAlgJ(Xs, theta, params, pcPWMp)
    % alpha(N, i, j) P(y_s_j=i| x_s_1, ...x_s_j, startT, T, PWMs)
    % scale(N, i) = P(x_s_i| startT, T, PWMs)
    % m x L
    kronMN = kron(1:params.m, ones(1, params.N));
    % the k+1 index is for base modes, 1 to k are for sub modes
    alphaSub = zeros(params.N, params.m, params.L + params.J, params.k);
    alphaBase = zeros(params.N, params.m, params.L);
    alphaT = zeros(params.N, params.m, params.L);
    scale = zeros(params.N, params.L);
    startE = matUtils.sumDim(theta.E, 2 : params.order);

    alphaBase(:, :, params.J+1) = (repmat(theta.startT, [1, params.N]) .* startE(:, Xs(:, 1))).';
    % alphaT(:, :, params.J+1) = alpha(:, :, params.J+1) * theta.T;
    scale(:, 1) = sum(alphaBase(:, :, params.J+1) + sum(alphaSub(:, :, params.J+1,:), 4), 2);
    if params.all2allMode
        theta.M = bsxfun(@times, theta.M, theta.F);
        theta.E = bsxfun(@times, theta.E, 1-theta.F);
    end

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
        assert(not(any(isnan(alphaBase(:)))))
        assert(not(any(isnan(alphaSub(:)))))
        fprintf('Forward algorithm %.2f%%\r', 100 * t / params.L);
        % N x m
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);

        % N x m
        if params.all2allMode
            newAlphaBase = alphaT(:, :, t - 1) .* Ep;
        else
            Tf = theta.T;
            Tf(eye(m)==1) = Tf(eye(m)==1) .* (1-theta.F);
            newAlphaBase = (alphaBase(:, :, t - 1) * Tf + sum(alphaSub(:,:,t-1,:), 4)) .* Ep;
        end


        % N x m x k
        alphaSlice = alphaT(:, :, t - theta.lengths);
        % N x m x k
        newAlphasSub = BaumWelchPWM.PWMstep(alphaSlice, Ms, t - theta.lengths' + 1, pcPWMp, params.J);

        % N x 1
        scale(:, t) = sum(newAlphaBase + sum(alphaSub, 3), 2);
        alphaBase(:, :, t) = bsxfun(@times, newAlphaBase, 1 ./ scale(:, t));
        alphaSub(:, :, t) = bsxfun(@times, newAlphasSub, 1 ./ scale(:, t));

        if params.all2allMode
            alphaT(:, :, t) = (alphaBase(:, :, t) + sum(alphaSub(:, :, t+J, :), 4)) * theta.T;
        else
            alphaT(:, :, t) = bsxfun(alphaBase(:, :, t), diag(theta.T)' .* theta.F);
        end
    end
    fprintf('\n');

    alphaSub = alphaSub(:, :, params.J+1:end);
    alphaSub = cat(3, alphaSub, zeros(params.N, params.m, params.J));
end


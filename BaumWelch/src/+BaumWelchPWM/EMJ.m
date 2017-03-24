
function [startT, T, Y, E, F, likelihood, gamma] = EMJ(Xs, m, maxIter, tEpsilon, order, PWMs, lengths)
    % PWMs - k x n x J emission matrix of m Jaspar PWMs with length J
    %        true length of i'th PWM< J and given in lengths(i) if a PWM is
    %        shorter than j, it is aligned to the end of the 3rd dimension.
    % Xs - N x L emission variables
    % m - amount of possible states (y)
    % n - amount of possible emissions (x)
    % maxIter - maximal iterations allowed
    % tEpsilon - probability to switch states will not exceed this number (if tEpsilon = 0.01,
    %            and m = 3 then the probability to stay in a state will not be less than 0.98)
    % order - the HMM order of the E matrix
    % initial estimation parameters
    [N,L] = size(Xs);

    [k, n, J] = size(PWMs);
    PWMsRep = permute(repmat(PWMs, [1, 1, 1, N]), [4, 3, 2, 1]);
    epsilon = 10 ^ -4;
    bestLikelihood = -Inf;
    repeat = 1;

    Xs1H = cat(2, zeros(N, J, n), matUtils.mat23Dmat(Xs, n));
    pcPWMp = preComputePWMp(PWMsRep, Xs1H, N, k, L);

    indices = reshape(matUtils.getIndices1D(Xs, order), [L-order+1, N]).';
    maxIndices = max(indices(:));
    indicesHotMap = matUtils.mat23Dmat(indices, maxIndices);
    indicesHotMap = cat(2, false(N, order-1, maxIndices), indicesHotMap);
    for rep = 1:repeat
        [startT, T, E, Y, F] = BaumWelchPWM.genRandParamsJ(m, n, order, k);

        iterLike = [];
        for it = 1:maxIter;
            % N x m x L
            % N x L
            % [alpha, scale] = BaumWelchPWM.forwardAlgJ(Xs, startT, T, Y, F, E, lengths, pcPWMp, J);
            % alpha = rand(N,m,L);
            % scale = rand(N,L);
            % N x m x L
            beta = BaumWelchPWM.backwardAlgJ(Xs, T, Y, F, E, scale, lengths, pcPWMp, J);

            % N x m x L
            % gamma_t(i) = P(y_t = i|x_1:L)
            gamma = alpha .* beta ./ repmat(sum(alpha .* beta, 2), [1, m, 1]);

            % N x m x L -> m x N x L
            [T, Y, F] = updateTYF(E, T, Y, F, Xs, alpha, beta, L, N, m, order);
            E = updateE(gamma, E, maxIndices, indicesHotMap);
            startT = updateStartT(gamma, startT);




            % T bound trick
            T = Tbound(T, tEpsilon, m);

            iterLike(end + 1) = sum(log(scale(:)));
            fprintf('Likelihood last iteration is %.2f', iterLike(end));
            if length(iterLike) > 1 && abs((iterLike(end) - iterLike(end-1)) / iterLike(end)) < epsilon
                % likelihood converged
                likelihood = iterLike(end);
                % fprintf('EM converged after %d iteratiosn: %f\n', it, likelihood);
                break
            end
        end % end of iterations loop

        if bestLikelihood < likelihood
            bestLikelihood = likelihood;
            bestE = E;
            bestT = T;
            bestGamma = gamma;
            bestStartT = startT;
        end
    end % end of repeat loop

    % return parameters with best likelihood
    E = bestE;
    T = bestT;
    gamma = bestGamma;
    startT = bestStartT;
    likelihood = bestLikelihood;
end

function newE = updateE(gamma, E, maxIndices, indicesHotMap)
    perGamma = permute(gamma, [2, 1, 3]);
    for i = 1:maxIndices
        % m x N x L -> m x 1
        E(:, i) = E(:, i) + sumDim(perGamma(:, indicesHotMap(:, :, i)), [2, 3]);
    end
    newE = bsxfun(@times, E, 1 ./ sum(E, order+1));
end
function newStartT = updateStartT(gamma, startT)
    startT = startT + mean(gamma(:, :, 1), 1).';
    % probability distribution normalization
    newStartT = startT / sum(startT);
end

function [newT, newY, newF] = updateTYF(E, T, Y, F, Xs, alpha, beta, L, N, m, order)
    TCorrection = zeros(m, m);
    YCorrection = zeros(m, k);
    FCorrection = zeros(m, 2);
    % todo: remove slow loop
    matSize = [m , 4 * ones(1, order)];
    kronMN = kron(1:m, ones(1, N));

    for t = 2 : L
        % for each letter in seqs we get more information about the transition matrix update
        % N x k
        PWMProb = getPWMp(BaumWelchPWM.PWMstep, Xs1H, t-1);
        % TODO: update y vector also
        % N x m x k
        M = repmat(alpha(:, :, t-1), [1,1,k]) .* beta(:, :, t + lengths - 1);
        % N x k x m
        M = bsxfun(@times, permute(M, [1,3,2]), PWMProb);
        % m x k x N
        M = bsxfun(@times, permute(M, [3,2,1]), Y);
        % m x k x N -> m x k
        % N x m
        Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % (m x N * N x m) .* (m x m)
        TCorrectiont = (alpha(:, :, t - 1).' * (beta(:, :, t) .* Ep)) .* T;
        TCorrection = TCorrection + (TCorrectiont / sum(sum(TCorrectiont)));

        YCorrectiont = sum(M, 3);
        YCorrection = YCorrection + (YCorrectiont / sum(sum(YCorrectiont)));

        FCorrectiont = [sumDim(M, [2, 3]), sum(TCorrectiont, 2)];
        FCorrection = FCorrection + FCorrectiont;
    end
    newT = T + TCorrection;
    newT = bsxfun(@times, newT, 1 ./ sum(newT, 2));
    newY = Y + YCorrection;
    newY = bsxfun(@times, newY, 1 ./ sum(newY, 2));
    newF = F + FCorrection;
    newF = bsxfun(@times, newF, 1 ./ sum(newF, 2));
    newF = newF(:, 1);
end

function newT = Tbound(T, tEpsilon, m)
    for i = 1 : m
        for j = 1 : m
            if i ~= j && T(i,j) > tEpsilon
                T(i, i) = T(i, i) + (T(i, j) - tEpsilon);
                T(i, j) = tEpsilon;
            end
        end
    end
    newT = T;
end


function pcPWMp = preComputePWMp(PWMsRep, Xs1H, N, k, L)
    % pcPWM - N x k x L-1
    PC_PWM_PROBABILITY_FILE = 'pcPWMp.mat';
    persistent p_pcPWMp
    if isempty(p_pcPWMp)
        try
            fprintf('Loading pre-computed PWM for sequences from %s...\n', PC_PWM_PROBABILITY_FILE);
            load(PC_PWM_PROBABILITY_FILE, 'pcPWMp')
            fprintf('Done.\n');
            p_pcPWMp = pcPWMp;
        catch
            fprintf('Pre-computing PWM probability on %d sequences', size(Xs1H, 1));
            pcPWMp = zeros(N, k, L-1);
            for t = 2:L
                fprintf('.');
                pcPWMp(:, :, t - 1) = BaumWelchPWM.getPWMp(PWMsRep, Xs1H, t);
            end
            fprintf('\n');
            save(PC_PWM_PROBABILITY_FILE, 'pcPWMp', '-v7.3')
            p_pcPWMp = pcPWMp;
        end
    end
    pcPWMp = p_pcPWMp;
end
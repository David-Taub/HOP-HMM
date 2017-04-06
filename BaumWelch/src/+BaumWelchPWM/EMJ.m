
function [startT, T, Y, E, F, likelihood, gamma] = EMJ(Xs, m, maxIter, tEpsilon, order, PWMs, lengths, pcPWMp)
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
    [N, L] = size(Xs);
    [k, n, J] = size(PWMs);

    epsilon = 10 ^ -4;
    bestLikelihood = -Inf;
    repeat = 1;

    % N x L -order + 1
    indices = reshape(matUtils.getIndices1D(Xs, order), [L-order+1, N]).';
    maxEIndex = max(indices(:));
    % N x L -order + 1 x maxEIndex
    indicesHotMap = matUtils.mat23Dmat(indices, maxEIndex);
    % N x L  x maxEIndex
    indicesHotMap = cat(2, false(N, order-1+J, maxEIndex), indicesHotMap);
    for rep = 1:repeat
        [startT, T, E, Y, F] = BaumWelchPWM.genRandParamsJ(m, n, order, k);

        iterLike = [];
        for it = 1:maxIter;
            % N x m x L
            % N x L
            [alpha, scale] = BaumWelchPWM.forwardAlgJ(Xs, startT, T, Y, F, E, lengths, pcPWMp, J);

            assert(not(any(isnan(alpha))))
            assert(not(any(isnan(scale))))
            % N x m x L
            beta = BaumWelchPWM.backwardAlgJ(Xs, T, Y, F, E, scale, lengths, pcPWMp, J);
            % beta = rand(N,m,L);
            assert(not(any(isnan(beta))))
            alpha = cat(3, alpha, zeros(N, m, J));
            beta = cat(3, beta, zeros(N, m, J));
            % N x m x L + J
            % gamma_t(i) = P(y_t = i|x_1:L)
            gamma = alpha .* beta ./ repmat(sum(alpha .* beta, 2), [1, m, 1]);
            assert(not(any(isnan(gamma))))

            % N x m x L -> m x N x L
            fprintf('Updating theta\n');
            [T, Y, F] = updateTYF(E, T, Y, F, Xs, alpha, beta, lengths, pcPWMp);
            E = updateE(gamma, E, maxEIndex, indicesHotMap);
            startT = updateStartT(gamma, startT);

            % T bound trick
            T = Tbound(T, tEpsilon, m);
            assert(not(any(isnan(T))))
            assert(not(any(isnan(Y))))
            assert(not(any(isnan(F))))
            assert(not(any(isnan(E))))
            assert(not(any(isnan(startT))))

            iterLike(end + 1) = sum(log(scale(:)));
            fprintf('Likelihood last iteration is %.2f\n', iterLike(end));
            if length(iterLike) > 1 && abs((iterLike(end) - iterLike(end-1)) / iterLike(end)) < epsilon
                % likelihood converged
                likelihood = iterLike(end);
                % fprintf('EM converged after %d iterations: %f\n', it, likelihood);
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

% indicesHotMap - N x L x maxEIndex
% gamma - N x m x L + J
% E - m x 4 x 4 x 4 x ... x 4 (order times)
function newE = updateE(gamma, E, maxEIndex, indicesHotMap)
    %  m x N x L + J
    fprintf('Update E\n');
    order = getOrder(E);
    perGamma = permute(gamma, [2, 1, 3]);
    for i = 1:maxEIndex
        % m x N x L -> m x 1
        E(:, i) = E(:, i) + sum(perGamma(:, indicesHotMap(:, :, i)), 2);
    end

    newE = bsxfun(@times, E, 1 ./ sum(E, order+1));
end

% gamma - N x m x L + J
function newStartT = updateStartT(gamma, startT)
    startT = startT + mean(gamma(:, :, 1), 1).';
    % probability distribution normalization
    newStartT = startT / sum(startT);
end

function order = getOrder(E)
    order = length(size(E)) - 1;
end

% Y - m x k
% F - k x 1
% T - m x m
% E - m x 4 x 4 x 4 x ... x 4 (order times)
% alpha - N x m x L + J
% beta - N x m x L + J
% lengths - k x 1
function [newT, newY, newF] = updateTYF(E, T, Y, F, Xs, alpha, beta, lengths, pcPWMp)
    order = getOrder(E);
    [m, k] = size(Y);
    [N, L] = size(Xs);
    TCorrection = zeros(m, m);
    YCorrection = zeros(m, k);
    FCorrection = zeros(m, 2);
    % todo: remove slow loop
    matSize = [m , 4 * ones(1, order)];
    kronMN = kron(1:m, ones(1, N));

    for t = 2 : L
        fprintf('\rUpdate TYF %.2f%%', 100*t/L);
        % for each letter in seqs we get more information about the transition matrix update
        % N x k
        PWMProb = pcPWMp(:, :, t);
        % N x m x k
        M = repmat(alpha(:, :, t-1), [1,1,k]) .* beta(:, :, t + lengths - 1);
        % N x k x m
        M = bsxfun(@times, permute(M, [1,3,2]), PWMProb);

        perY = permute(Y, [3, 2, 1]);
        for i = 1:N
            M(i, :, :) = M(i, :, :) .* perY;
        end

        % todo,  delete:
        % % N x k x m -> m x k x N
        % M = permute(M, [3,2,1]);
        % % Y - m x k
        % M = bsxfun(@times, M, Y);
        % YCorrectiont = sum(M, 3);

        % m x k x N -> m x k
        % N x m
        Ep = BaumWelchPWM.getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % (m x N * N x m) .* (m x m)
        TCorrectiont = (alpha(:, :, t - 1)' * (beta(:, :, t) .* Ep)) .* T;
        TCorrection = TCorrection + (TCorrectiont / sum(sum(TCorrectiont)));


        YCorrectiont = matUtils.sumDim(M, 1)';
        YCorrection = YCorrection + (YCorrectiont / sum(sum(YCorrectiont)));

        FCorrectiont = [matUtils.sumDim(M, [1, 2]), sum(TCorrectiont, 2)];
        FCorrection = FCorrection + FCorrectiont;
    end
    fprintf('\n');
    newT = T + TCorrection;
    newT = bsxfun(@times, newT, 1 ./ sum(newT, 2));
    newY = Y + YCorrection;
    newY = bsxfun(@times, newY, 1 ./ sum(newY, 2));
    newF = [F, 1-F] + FCorrection;
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


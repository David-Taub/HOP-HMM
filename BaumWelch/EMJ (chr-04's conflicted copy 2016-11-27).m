    
function [startT, T, Y, E, likelihood, gamma] = EMJ(Xs, m, maxIter, tEpsilon, order, PWMs, lengths)
    % PWMs - k x n x J emission matrix of m Jaspar PWMs with length J 
    %        true length of i'th PWM< J and given in lengths(i) if a PWM is 
    %        shorter than j, it is aligned to the end of the 3rd dimension.
    % Xs - N x L emission variables
    % m - ammount of possible states (y)
    % n - amount of possible emmissions (x)
    % maxIter - maximal iteraions allowed
    % tEpsilon - probability to switch states will not exceed this number (if tEpsilon = 0.01, 
    %            and m = 3 then the probability to stay in a state will not be less than 0.98)
    % order - the HMM order of the E matrix
    % initial estimation parameters
    [N,L] = size(Xs);

    [k, n, J] = size(PWMs);
    Xs1HRep = repmat(mat23Dmat(Xs, n), [1, 1, 1, k]);
    PWMsRep = permute(repmat(PWMs, [1, 1, 1, N]), [4, 3, 2, 1]);
    subModePost = zeros(N, k, L-J+1);
    for t = 1:L-J+1
        % N x k
        subModePost(:, :, t) = sumDim(PWMsRep .* Xs1HRep(:, t:t+J-1, :, :), [2,3]);
    end
    epsilon = 10 ^ -4;
    bestLikelihood = -Inf;
    repeat = 1;
    indices = reshape(getIndeices1D(Xs, order), [L-order+1, N]).';
    maxIndices = max(indices(:)); %might be lower than n^order therefore shorter loop
    indicesHotMap = mat23Dmat(indices, maxIndices);
    % N x L x maxIndices
    indicesHotMap = cat(2, false(N, order-1, maxIndices), indicesHotMap);
    for rep = 1:repeat
        [startT, T, E, Y] = genRandParamsJ(m, n, order, k);

        iterLike = [];
        for it = 1:maxIter;
            % N x m x L
            % m x L
            submodeLike
            [alpha, scale] = forwardAlgJ(Xs, startT, T, Y, PWMsRep, lengths);
            % N x m x L
            beta = backwardAlgJ(Xs, startT, T, E, Y, scale, PWMsRep, lengths);

            % N x m x L
            % gamma_t(i) = P(y_t = i|x_1:L)
            gamma = alpha .* beta ./ repmat(sum(alpha .* beta, 2), [1, m, 1]);
            
            % update E and T
            T = updateT(E, T, Xs, alpha, beta, L, N, m, order);
            
            % N x m x L -> m x N x L
            perGamma = permute(gamma, [2, 1, 3]);
            %todo: consider iteration over L, might be faster if order > 4~5
            for i = 1:maxIndices
                % N x L x maxIndices
                % m x N x L -> m x 1
                EUpdate = sumDim(perGamma(:, indicesHotMap(:, :, i)), [2, 3]);
                E(:, i) = E(:, i) + EUpdate;
            end
            startT = startT + mean(gamma(:, :, 1), 1).';

            % probability distribution normalization
            startT = startT / sum(startT);
            T = bsxfun(@times, T, 1 ./ sum(T, 2));
            E = bsxfun(@times, E, 1 ./ sum(E, order+1));


            % T bound trick
            T = Tbound(T, tEpsilon, m);
            iterLike(end + 1) = sum(log(scale(:)));
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
% TODO: get PWMs also
function newT = updateT(E, T, Xs, alpha, beta, L, N, m, order)
    xi = zeros(m, m);
    % todo: remove slow loop
    matSize = [m , 4 * ones(1, order)];
    kronMN = kron(1:m, ones(1, N));
    for t = 2 : L
        % N x m
        Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % for each letter in seqs we get more information about the transition matrix update
        % (m x N * N x m) .* (m x m)
        newXi = (alpha(:, :, t-1).' * (beta(:, :, t) .* Ep)) .* T;
        xi = xi + (newXi / sum(sum(newXi)));
    end
    newT = T + xi;
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
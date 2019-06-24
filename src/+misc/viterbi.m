% X - N x L
% Y - N x L x 2
% pcPWMp - N x k x L
function Y = viterbi(params, theta, X, pcPWMp)
    [N, L] = size(X);
    fprintf('Running Viterbi on %d sequences of length %d\n', N, L);
    pcPWMp = cat(3, -inf(N, params.k, params.J), pcPWMp);
    Y = zeros(N, L, 2);
    matSize = [params.m , params.n * ones(1, params.order)];
    kronMN = kron(1:params.m, ones(1, N));
    % N x m x L
    Eps = EM.getEp3d(theta, params, X, 1:L);
    for i = 1:N
        O1 = -inf(params.m, L+params.J);
        O2 = zeros(params.m, L, 2);
        % m x L x 2
        O1(:, 1+params.J) = theta.startT + Eps(i, :, 1)';
        % Viterbi dynamic programming
        for t = 2:L
            % m x m
            baseSteps = repmat(O1(:, t-1+params.J), [1, params.m]) + theta.T + repmat(Eps(i, :, t), [params.m, 1]);
            % m x k
            subSteps = zeros(params.m, params.k);
            for j = 1:params.k
                subSteps(:, j) = O1(:, t-params.lengths(j)-1+params.J) + permute(pcPWMp(i, j, t-params.lengths(j)'+params.J), [1,3,2]);
            end
            % m x k
            subSteps = subSteps + theta.G + repmat(Eps(i, :, t)', [1, params.k]);
            O2(:, t, 1) = 1:params.m;
            [O1(:, t+params.J), O2(:, t, 2)] = max(subSteps, [], 2);
            % 1 x m
            [maxBaseStep, maxBaseStepInds] = max(baseSteps, [], 1);
            baseMask = O1(:, t+params.J) <= maxBaseStep';

            % if baseMask(1) == 0
            %     keyboard
            % end
            O2(baseMask, t, 2) = 0;
            O2(baseMask, t, 1) = maxBaseStepInds(baseMask);
            O1(baseMask, t+params.J) = maxBaseStep(baseMask);
        end

        % back trailing
        Y(i, L, 2) = 0;
        [~, Y(i, L, 1)] = max(O1(:, end), [], 1);
        t = L;
        while t >= 2
            % m x m
            Y(i, t-1, :) = O2(Y(i, t, 1), t, :);
            if Y(i, t-1, 2) > 0
                motifLen = params.lengths(Y(i, t-1, 2));
                Y(i, t-motifLen-1:t-1, :) = repmat(Y(i, t-1, :), [1, motifLen+1, 1]);
                Y(i, t-motifLen-1, 2) = 0;
                t = t-params.lengths(Y(i, t-1, 2))-1;
            else
                t = t - 1;
            end
        end
    end
end
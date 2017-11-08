% X - N x L
% Y - N x L
% pcPWMp - N x k x L
function Y = viterbi(theta, params, X, pcPWMp)
    [N, L] = size(X);
    % N x k x J + L
    pcPWMp = cat(3, -inf(N, params.k, params.J), pcPWMp);
    Y = zeros(N, L, 2);
    % N x m x L
    Eps = EM.getEp3d(theta, params, X, 1:L);
    tt = []; tt2 = [];
    for i = 1:N
        fprintf('%d / %d\r',i, N)
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
                subSteps(:, j) = pcPWMp(i, j, t-params.lengths(j)+params.J);
            end
            % m x k
            subSteps = subSteps + O1(:, t-params.lengths'-1+params.J) + theta.G + repmat(Eps(i, :, t)', [1, params.k]);

            O2(:, t, 1) = 1:params.m;
            [O1(:, t+params.J), O2(:, t, 2)] = max(subSteps, [], 2);
            % 1 x m
            [maxBaseStep, maxBaseStepInds] = max(baseSteps, [], 1);
            O1(:, t+params.J)
            maxBaseStep
            tt = [tt;maxBaseStep];
            tt2 = [tt2;O1(:, t+params.J)'];
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
        while t > 2
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

        subplot(2,3,1); imagesc(O1); colorbar;
        subplot(2,3,2); imagesc(O2(:, :, 1)); colorbar;
        subplot(2,3,3); imagesc(O2(:, :, 2)); colorbar;
        subplot(2,3,4); imagesc(Y(1:i, :, 1)); colorbar;
        subplot(2,3,5); imagesc(Y(1:i, :, 2)); colorbar;
        subplot(2,3,6); plot(tt(:, 2)'); hold on; plot(tt2(:, 2)'); hold off;

        drawnow
        keyboard
    end
end
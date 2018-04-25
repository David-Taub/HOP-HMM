
function [bestTheta, bestLikelihood] = EM(dataset, params, maxIter, doResample, doESharing)
    % X - N x L emission variables
    % m - amount of possible states (y)
    % n - amount of possible emissions (x)
    % maxIter - maximal iterations allowed
    % tEpsilon - probability to switch states will not exceed this number (if tEpsilon = 0.01,
    %            and m = 3 then the probability to stay in a state will not be less than 0.98)
    % pcPWMp - N x k x L
    % order - the HMM order of the E matrix
    % initial estimation parameters

    [N, L] = size(dataset.X);
    fprintf('Starting EM algorithm on %d x %d\n', N, L)
    bestLikelihood = -Inf;
    repeat = 1;
    % N x L - order + 1
    indices = reshape(matUtils.getIndices1D(dataset.X, params.order, params.n), [L-params.order+1, N]).';
    % N x L - order + 1 x maxEIndex
    indicesHotMap = matUtils.mat23Dmat(indices, params.n ^ params.order);
    % N x L  x maxEIndex
    indicesHotMap = cat(2, false(N, params.order-1, params.n ^ params.order), indicesHotMap);
    % figure
    for rep = 1:repeat
        % X = X(randperm(N), :);
        initTheta = misc.genTheta(params);
        [iterLike, theta] = singleRunEM(dataset, params, initTheta, maxIter, indicesHotMap, N, L, doResample, doESharing);
        if bestLikelihood < iterLike(end)
            bestLikelihood = iterLike(end);
            bestTheta = theta;
        end
    end
end

function [iterLike, theta] = singleRunEM(dataset, params, initTheta, maxIter, indicesHotMap, N, L, doResample, doESharing)
    LIKELIHOOD_THRESHOLD = 10 ^ -4;
    theta = initTheta;
    iterLike = [];
    for it = 1:maxIter
        tic
        % alpha - N x m x L
        % beta - N x m x L
        % pX - N x 1
        % xi - N x m x m x L
        % gamma - N x m x L
        % psi - N x m x k x L
        [alpha, beta, pX, xi, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
        close all;
        % show.showTheta(theta);
        % show.seqSampleCertainty(params, dataset.Y, gamma, psi, 8, true);
        theta.E = updateE(gamma, params, indicesHotMap);
        % drawStatus(theta, params, alpha, beta, gamma, pX, xi, psi);
        if doESharing
            theta.E = log(repmat(mean(exp(theta.E), 1), [params.n, ones(1, params.order)]));
        end
        % fprintf('Update G\n');
        fprintf('. ');
        [theta.G, theta.T] = updateGT(params, theta, xi, gamma, psi, doResample);
        % fprintf('Update startT\n');
        theta.startT = updateStartT(gamma);
        iterLike(end+1) = matUtils.logMatSum(pX, 1);% / (N*L);
        assert(not(any(isnan(theta.T(:)))))
        assert(not(any(isnan(theta.E(:)))))
        assert(not(any(isnan(theta.G(:)))))
        assert(not(any(isnan(alpha(:)))))
        assert(not(any(isnan(beta(:)))))
        motifsPer = sum(exp(theta.G(:)), 1).*100;
        timeLapse = toc();
        R = cov(theta.G);
        fprintf('It %d: log-like: %.2f Time: %.2fs, motifs: ~%.2f%% cov: %.2f\n', it, iterLike(end), timeLapse, motifsPer, mean(R(:)));
        if length(iterLike) > 1 && abs((iterLike(end) - iterLike(end-1)) / iterLike(end)) < LIKELIHOOD_THRESHOLD
            fprintf('Converged\n');
            break
        end
    end % end of iterations loop
 end


% psi - N x m x k x L
% pX - N x 1
% gamma - N x m x L
function drawStatus(theta, params, alpha, beta, gamma, pX, xi, psi)
    figure
    YsEst = mean(gamma, 3);
    YsEst = YsEst(:, 1);
    subplot(2,2,1);
    scatter(1:size(pX, 1), YsEst); colorbar;
    title('Px')
    % [~, YsEst] = max(alpha, [], 2);
    % subplot(2,3,2);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % % title('alpha')
    % [~, YsEst] = max(beta(:,:,1:end), [], 2);
    % subplot(2,3,3);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % title('beta')
    subplot(2,2,2);plot(exp(theta.E(:)));
    plot(exp(theta.E(:)));
    ylim([0,1]);
    title('E')

    subplot(2,2,3);
    plot(permute(gamma(1,1,:), [3,2,1]));
    hold on;
    % N x m x k x L
    plot(permute(matUtils.logMatSum(psi(1,1,:,:), 3), [4,3,2,1]));
    title('Postirior')
    subplot(2,2,4);plot(exp(theta.G(:)));
    ylim([0,1]);
    title('G')
    drawnow
end

% indicesHotMap - N x L x maxEIndex
% gamma - N x m x L
% E - m x 4 x 4 x 4 x ... x 4 (order times)
function newE = updateE(gamma, params, indicesHotMap)
    %  m x N x L + J
    % m x N x L
    perGamma = permute(gamma, [2, 1, 3]);
    newE = -inf([params.m, params.n * ones(1, params.order)]);
    for i = 1:(params.n ^ params.order)
        % m x N x L -> m x 1
        newE(:, i) = matUtils.logMatSum(perGamma(:, indicesHotMap(:, :, i)), 2);
    end
    newE = matUtils.logMakeDistribution(newE);
end


% gamma - N x m x L
function newStartT = updateStartT(gamma)
    newStartT = matUtils.logMakeDistribution(matUtils.logMatSum(gamma(:, :, 1), 1))';
    % probability distribution normalization
    assert(size(newStartT, 1) > 0)
end



% alpha - N x m x L
% beta - N x m x L
% gamma - N x m x L
% xi - N x m x m x L
% psi - N x m x k x L
% newG - m x k
% newT - m x m
function [newG, newT] = updateGT(params, theta, xi, gamma, psi, doResample)
    [N, ~, L] = size(gamma);

    % note: batch trick is used to reduce the
    % calculation errors due to summing
    % many very small numbers in log space

    % N x m x k x L
    batchAmount = ceil(N / params.batchSize);
    % N x m x k+m x L
    psiXiMerged = cat(3, psi, xi);
    psiXiMerged = psiXiMerged(:, :, :, 1:end-params.J);
    % psiXiMerged = psiXiMerged - permute(repmat(gamma(:,:,1:end-params.J), [1, 1, 1, params.k+params.m]), [1,2,4,3]);

    psiXiMerged = matUtils.logMatSum(psiXiMerged, 4);
    mergedBatches = -inf(batchAmount, params.m, params.k+params.m);
    for u = 1:batchAmount
        batchMask = mod(1:N, batchAmount) == u-1;
        batch = psiXiMerged(batchMask, :, :);
        % batchSize x m x k+m x L-J
        mergedBatches(u, :, :) = matUtils.logMatSum(batch, 1);
        mergedBatches(u, :, :) = matUtils.logMakeDistribution(mergedBatches(u, :, :));
    end

    % batchSize x m x k+m -> m x k+m
    mergedAveraged = matUtils.logMatSum(mergedBatches, 1);
    mergedAveraged = matUtils.logMakeDistribution(mergedAveraged);
    mergedAveraged = permute(mergedAveraged, [2,3,1]);
    G = mergedAveraged(:, 1:params.k);
    T = mergedAveraged(:, params.k+1:end);
    [newG, newT] = EM.GTbound(params, G, T, doResample);
end


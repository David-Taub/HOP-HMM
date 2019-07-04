function [theta, iterLike] = EMIteration(params, dataset, inputTheta, doGTBound)
    N = size(dataset.X, 1);
    LEARNING_RATE = 0.7;
    if ~isfield(dataset, 'XIndicesHotMap')
        % N x L - order + 1
        dataset.XIndicesHotMap = misc.genXInidcesHotMap(params, dataset);
    end
    batchAmount = floor(N / params.batchSize);
    batchesTheta = inputTheta;
    batchesTheta.T = batchesTheta.T * 0;
    batchesTheta.E = batchesTheta.E * 0;
    batchesTheta.G = batchesTheta.G * 0;
    batchesTheta.startT = batchesTheta.startT * 0;
    batchesLikelihood = 0;
    for u = 1:batchAmount
        % N x m x k+m x L
        batchMask = mod(1:N, batchAmount) == u - 1;
        fprintf('Batch %d / %d [%d / %d] \n', u, batchAmount, sum(batchMask, 2), N);
        XBatch = dataset.X(batchMask, :);
        pcPWMpBatch = dataset.pcPWMp(batchMask, :, :);
        XIndicesHotMapBatch = dataset.XIndicesHotMap(batchMask, :, :);

        % alpha - N x m x L
        % beta - N x m x L
        % pX - N x 1
        % xi - N x m x m x L
        % gamma - N x m x L
        % psi - N x m x k x L
        [alpha, beta, pX, xi, gamma, psi] = EM.EStep(params, inputTheta, XBatch, pcPWMpBatch);
        % close all;
        % show.showTheta(theta);
        % show.seqSampleCertainty(params, dataset.Y, gamma, psi, 8, true);
        E = updateE(gamma, params, XIndicesHotMapBatch);
        % drawStatus(theta, params, alpha, beta, gamma, pX, xi, psi);
        if params.doESharing
            % theta.E = log(repmat(mean(exp(theta.E), 1), [params.n, ones(1, params.order)]));
            E(1:end - params.backgroundAmount, :) = repmat(log(mean(exp(E(1:end - params.backgroundAmount, :)), 1)), [params.n - params.backgroundAmount, 1]);
        end
        fprintf('. ');
        startT = updateStartT(gamma);
        [G, T] = updateGT(params, xi, gamma, psi);
        assert(not(any(isnan(startT(:)))));
        assert(not(any(isnan(T(:)))));
        assert(not(any(isnan(E(:)))));
        assert(not(any(isnan(G(:)))));
        assert(not(any(isnan(alpha(:)))));
        assert(not(any(isnan(beta(:)))));
        batchesTheta.T = batchesTheta.T + exp(T);
        batchesTheta.E = batchesTheta.E + exp(E);
        batchesTheta.G = batchesTheta.G + exp(G);
        batchesTheta.startT = batchesTheta.startT + exp(startT);
        batchesLikelihood = batchesLikelihood + sum(pX, 1);
    end
    % geometric average is more stable than regular mean for very small likelihood values
    iterLike = batchesLikelihood / N;
    theta.T = log(exp(inputTheta.T) .* (1 - LEARNING_RATE) + LEARNING_RATE .* batchesTheta.T / batchAmount);
    theta.E = log(exp(inputTheta.E) .* (1 - LEARNING_RATE) + LEARNING_RATE .* batchesTheta.E / batchAmount);
    theta.G = log(exp(inputTheta.G) .* (1 - LEARNING_RATE) + LEARNING_RATE .* batchesTheta.G / batchAmount);
    theta.startT = log(exp(inputTheta.startT) .* (1 - LEARNING_RATE) + LEARNING_RATE .* batchesTheta.startT / batchAmount);
    if doGTBound
        [theta.G, theta.T] = EM.GTbound(params, theta.G, theta.T);
    end
    if params.doResampling
        [theta.E, theta.G] = EM.resampleEG(params, theta.E, theta.G);
    end
    % clf
    % show.showTwoThetas(params, dataset.theta, theta, false, sprintf('%d', it), 'tmp.jpg');
    % drawnow;
end


% dataset.XIndicesHotMap - N x L x maxEIndex
% gamma - N x m x L
% E - m x 4 x 4 x 4 x ... x 4 (order times)
function newE = updateE(gamma, params, XIndicesHotMap)
    %  m x N x L + J
    % m x N x L
    permGamma = permute(gamma, [2, 1, 3]);
    newE = -inf([params.m, params.n * ones(1, params.order)]);
    for i = 1:(params.n ^ params.order)

        % m x N x L -> m x 1
        kmerPosteriorProb = permGamma(:, XIndicesHotMap(:, :, i));
        if size(kmerPosteriorProb, 2) > 0
            newE(:, i) = matUtils.logMatSum(kmerPosteriorProb, 2);
        end
    end
    newE = matUtils.logMakeDistribution(newE);
end


% gamma - N x m x L
function newStartT = updateStartT(gamma)
    newStartT = matUtils.logMakeDistribution(matUtils.logMatSum(gamma(:, :, 1), 1))';
    % probability distribution normalization
    assert(size(newStartT, 1) > 0);
end


% alpha - N x m x L
% beta - N x m x L
% gamma - N x m x L
% xi - N x m x m x L
% psi - N x m x k x L
% newG - m x k
% newT - m x m
function [newG, newT] = updateGT(params, xi, gamma, psi)
    [N, ~, L] = size(gamma);

    % N x m x k+m x L
    psiXiMerged = cat(3, psi, xi);
    psiXiMerged = psiXiMerged(:, :, :, 1:end-params.J);
    % psiXiMerged = psiXiMerged - permute(repmat(gamma(:,:,1:end-params.J), [1, 1, 1, params.k+params.m]), [1,2,4,3]);

    psiXiMerged = matUtils.logMatSum(psiXiMerged, 4);
    mergedAveraged = matUtils.logMatSum(psiXiMerged, 1);
    mergedAveraged = matUtils.logMakeDistribution(mergedAveraged);
    mergedAveraged = permute(mergedAveraged, [2,3,1]);
    newG = mergedAveraged(:, 1:params.k);
    newT = mergedAveraged(:, params.k+1:end);
end

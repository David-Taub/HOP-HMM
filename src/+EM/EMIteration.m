function [theta, iterLike] = EMIteration(params, dataset, inputTheta, doGTBound)
    [N, L] = size(dataset.X);
    % the bigger this number is, the more time a batch takes
    % EXPECTED_COMPUTATION_IN_BATCH = 50000;
    % EXPECTED_COMPUTATION_IN_BATCH = 100000;
    % EXPECTED_COMPUTATION_IN_BATCH = 150000; % 20-40k vars a second
    % EXPECTED_COMPUTATION_IN_BATCH = 200000;
    EXPECTED_COMPUTATION_IN_BATCH = 500000; % 23k vars a second
    % EXPECTED_COMPUTATION_IN_BATCH = 1000000;
    batchSize = min(N, ceil(EXPECTED_COMPUTATION_IN_BATCH / L));
    batchAmount = floor(N / batchSize);
    assert(batchAmount > 0)
    batchesTheta = inputTheta;
    batchesTheta.T = zeros(params.m, params.m);
    batchesTheta.E = zeros([params.m, ones(1, params.order) * params.n]);
    batchesTheta.G = zeros(params.m, params.k);;
    batchesTheta.startT = zeros(params.m, 1);;
    batchesLikelihood = 0;
    for currentBatchIndex = 1:batchAmount
        t = tic();
        % N x m x k+m x L
        batchMask = mod(1:N, batchAmount) == currentBatchIndex - 1;
        fprintf('Batch %d / %d [%d / %d] ', currentBatchIndex, batchAmount, sum(batchMask, 2), N);
        XBatch = dataset.X(batchMask, :);
        pcPWMpBatch = dataset.pcPWMp(batchMask, :, :);
                % N x L - order + 1
        XIndicesHotMapBatch = misc.genXInidcesHotMap(params, dataset.X(batchMask, :));

        % alpha - N x m x L
        % beta - N x m x L
        % pX - N x 1
        % xi - N x m x m x L
        % gamma - N x m x L
        % psi - N x m x k x L
        [alpha, beta, pX, xi, gamma, psi] = EM.EStep(params, inputTheta, XBatch, pcPWMpBatch);
        E = updateE(gamma, params, XIndicesHotMapBatch);
        assert(not(any(isnan(E(:)))));
        if params.doESharing
            % theta.E = log(repmat(mean(exp(theta.E), 1), [params.n, ones(1, params.order)]));
            E(:, :) = repmat(log(mean(exp(E(:, :)), 1)), [params.m, 1]);
        end
        fprintf('. ');
        startT = updateStartT(gamma);
        [G, T] = updateGT(params, xi, gamma, psi);
        batchesTheta.T = batchesTheta.T + exp(T);
        batchesTheta.E = batchesTheta.E + exp(E);
        batchesTheta.G = batchesTheta.G + exp(G);
        batchesTheta.startT = batchesTheta.startT + exp(startT);
        batchesLikelihood = batchesLikelihood + sum(pX, 1);
        assert(not(any(isnan(batchesTheta.startT(:)))));
        assert(not(any(isnan(batchesTheta.T(:)))));
        assert(not(any(isnan(batchesTheta.E(:)))));
        assert(not(any(isnan(batchesTheta.G(:)))));
        batchTime = toc(t);
        fprintf('[%.2f secs, %.2f vars/sec]\n', batchTime, L * batchSize / batchTime);
    end
    % geometric average is more stable than regular mean for very small likelihood values
    iterLike = batchesLikelihood / N;
    theta.T = log(exp(inputTheta.T) .* (1 - params.learningRate) + params.learningRate .* batchesTheta.T / batchAmount);
    theta.E = log(exp(inputTheta.E) .* (1 - params.learningRate) + params.learningRate .* batchesTheta.E / batchAmount);
    theta.G = log(exp(inputTheta.G) .* (1 - params.learningRate) + params.learningRate .* batchesTheta.G / batchAmount);
    theta.startT = log(exp(inputTheta.startT) .* (1 - params.learningRate) + params.learningRate .* batchesTheta.startT / batchAmount);
    % theta.startT = inputTheta.startT;
    if doGTBound > 0
        theta = EM.GTbound3(params, theta);
    end
    if params.doResampling
        [theta.E, theta.G, theta.T] = EM.resampleEG(params, theta.E, theta.G, theta.T);
    end
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

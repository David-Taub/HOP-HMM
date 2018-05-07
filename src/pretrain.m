
% doResample - if True will resample G of mode if is too similar to another mode (uses threashold)
% doESharing - Each EM iteration, averaging the E across all modes, and using the average in all modes

function pretrain(mergedPeaksMin, k, tissueList, backgroundIndex)
    dbstop if error
    close all;
    doResample = false;
    doESharing = false;
    m = 2;
    params = misc.genParams(m, k);
    params.NperTissue = 1000;
    testTrainRatio = 0.15;
    paramsTotal = params
    paramsTotal.m = length(tissueList);

    totalTheta = misc.genTheta(paramsTotal);
    for i = 1:length(tissueList)
        dataset = preprocess(params, mergedPeaksMin, tissueList(i), backgroundIndex);
        [theta, ~] = EM.EM(dataset, params, 50, doResample, doESharing);
        totalTheta.E(i, :) = theta.E(1, :);
        totalTheta.G(i, :) = theta.G(1, :);
        show.showTheta(totalTheta);
        [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
        show.seqSampleCertainty(params, dataset.Y, gamma, psi, 8, false);
        classify(theta, params, dataset);
    end
    E = totalTheta.E;
    G = totalTheta.G;
    save('../data/precomputation/pretrainedTheta.mat', 'tissueList', 'E', 'G');
end


function [dataset] = preprocess(params, mergedPeaksMin, tissueIndex, backgroundIndex)
    L = size(mergedPeaksMin.seqs, 2);
    X = mergedPeaksMin.seqs;

    if isfield(mergedPeaksMin, 'overlaps')
        overlaps = mergedPeaksMin.overlaps;
        mask = sum(overlaps > 0, 2) == 1;
        mask1 = mask & overlaps(:, tissueIndex) > 0;
        mask2 = mask & overlaps(:, backgroundIndex) > 0;
        pN = min(sum(mask1, 1), sum(mask2, 1));
        inds = find(mask1);
        mask1 = false(size(overlaps, 1), 1);
        mask1(inds(1:pN)) = true;
        inds = find(mask2);
        mask2 = false(size(overlaps, 1), 1);
        mask2(inds(1:pN)) = true;

        mask = mask & (mask1 | mask2);
        X = X(mask, :);
        overlaps = overlaps(mask, :);
        overlaps = overlaps(:, sort([tissueIndex, backgroundIndex]));
        Y = repmat((overlaps > 0) * [1;2], [1, L]);
        assert(not(any(Y(:)==0)))
    end
    % N x k x L
    dataset.title = 'dataset';
    dataset.X = X;
    dataset.Y = Y;
    dataset.pcPWMp = misc.preComputePWMp(X, params);
    % if isfield(mergedPeaksMin, 'Y')
    %     dataset.Y = mergedPeaksMin.Y(datasetMask, :);
    % end
    % if isfield(mergedPeaksMin, 'Y2')
    %     mergedPeaksMin.Y2
    %     dataset.Y2 = mergedPeaksMin.Y2(trainMask, :);
    % end
end

function perm = findCorrectPermute(params, Y, YEst)
    % N * L x m
    Y1Hot = matUtils.vec2mat(Y(:)', params.m)';
    YEst1Hot = matUtils.vec2mat(YEst(:)', params.m)';
    % 1 x m x N * L
    Y1HotPer = permute(Y1Hot, [3, 2, 1]);
    % m x 1 x N * L
    YEst1HotPer = permute(YEst1Hot, [2, 3, 1]);
    distMat = sum(repmat(Y1HotPer, [params.m, 1, 1]) == repmat(YEst1HotPer, [1, params.m, 1]), 3);
    perm = zeros(params.m, 1);
    for i = 1:params.m
        [~, I] = max(distMat(:));
        [I_row, I_col] = ind2sub(size(distMat),I);
        perm(I_col) = I_row;
        distMat(I_row, :) = -1;
        distMat(:, I_col) = -1;
    end
end


% Y - N x L
% YEst - N x L
function theta = permuteTheta(theta, params, Y, YEst)
    perm = findCorrectPermute(params, Y, YEst);
    theta.T = theta.T(perm, :);
    theta.T = theta.T(:, perm);
    theta.G = theta.G(perm, :);
    theta.E(:, :) = theta.E(perm, :);
    theta.startT = theta.startT(perm);
end

% P - u x 1
% Q - u x 1
function ret = relativeEntropy(P, Q)
    ret = sum(P .* (log(P) - log(Q)), 1);
end

% Y - N x L
function YEstViterbi = classify(theta, params, dataset)
    [N, L] = size(dataset.X);
    % N x m x L
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
    % EM.drawStatus(theta, params, gamma);
    YEstViterbi = misc.viterbi(theta, params, dataset.X, dataset.pcPWMp);
    YEstViterbiAcc = YEstViterbi(:, :, 1) == dataset.Y; YEstViterbiAcc = sum(YEstViterbiAcc(:)) ./ length(YEstViterbiAcc(:));

    YEstMax = maxPostEstimator(theta, params, psi, gamma);
    YEstMaxAcc = YEstMax(:, :, 1) == dataset.Y; YEstMaxAcc = sum(YEstMaxAcc(:)) ./ length(YEstMaxAcc(:));

    figure
    subplot(1,6,1);
    imagesc(YEstViterbi(:,:,1)'); colorbar; title(['Viterbi States', num2str(YEstViterbiAcc)]);
    subplot(1,6,2);
    imagesc(YEstViterbi(:,:,2)'); colorbar;title('Viterbi Motifs');
    subplot(1,6,3);
    imagesc(YEstMax(:,:,1)'); colorbar; title(['MaxPosterior States', num2str(YEstMaxAcc)]);
    subplot(1,6,4);
    imagesc(YEstMax(:,:,2)'); colorbar;title('MaxPosterior Motifs');
    subplot(1,6,5);
    imagesc(dataset.Y'); colorbar;title(['Real States (', dataset.title, ')']);
    if isfield(dataset, 'Y2')
        subplot(1,6,6);
        imagesc(dataset.Y2'); colorbar;title('Real Motifs');
    end

end

function YEst = maxPostEstimator(theta, params, psi, gamma)
    [N, ~, L] = size(gamma);
    % N x L x m
    gammaPer = permute(gamma, [1, 3, 2]);
    [gammaMaxVals, YEstBaseStates] = max(gammaPer, [], 3);

    % N x m x k x L -> N x L x m x k
    psiPer = cat(2, -inf(N, params.J, params.m, params.k), permute(psi, [1, 4, 2, 3]));
    skewedPsi = -inf(N, L, params.m, params.k);
    for l = 1:params.k
        for u = 1:params.lengths(l)
            skewedPsi(:, :, :, l) = max(skewedPsi(:, :, :, l), psiPer(:, [1:L] + params.J - u, :, l));
        end
    end
    % N x L x m x k -> N x L
    [psiMaxVals, psiMaxInd] = max(skewedPsi(:,:,:), [], 3);
    subStates = floor((psiMaxInd - 1) / params.m) + 1;
    baseStates = mod(psiMaxInd - 1, params.m) + 1;
    subStateMask = psiMaxVals > gammaMaxVals;
    % N x L
    YEstBaseStates(subStateMask) = baseStates(subStateMask);
    YEstSubStates = zeros(N, L);
    YEstSubStates(subStateMask) = subStates(subStateMask);
    YEst = cat(3, YEstBaseStates, YEstSubStates);
end


% posterior - N x m x L
function posterior = calcPosterior(params, gamma, psi)
    N = size(gamma, 1);
    posterior = gamma;
    for l = 1:params.k
        for t = 1:params.lengths(l)
            subStatePost = permute(cat(4, -inf(N, params.m, 1, t), psi(:, :, l, 1:end-t)), [1,2,4,3]);
            posterior = matUtils.logAdd(posterior, subStatePost);
        end
    end
    posterior = exp(posterior);
end


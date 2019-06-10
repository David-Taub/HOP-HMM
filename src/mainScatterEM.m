
function mainScatterEM()
    conf.m = 5
    conf.backgroundAmount = 0
    conf.k = 20
    conf.doResample = false
    conf.doESharing = false
    conf.doBound = true
    conf.N = 3000;
    conf.L = 1000;
    conf.maxIters = 1000;
    conf.canCrossLayer = true;
    conf.patience = 3;
    main(conf);
end

% doResample - if True will resample G of mode if is too similar to another mode (uses threshold)
% doESharing - Each EM iteration, averaging the E across all modes, and using the average in all modes
function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.canCrossLayer);
    [test, train] = crossValidationSplit(params, mergedPeaksMin, 0.01);
    [thetaEst, ~] = EM.EM(train, params, conf.maxIters, conf.doResample, conf.doESharing, conf.doBound, conf.patience);
    thetaOrig = mergedPeaksMin.theta;
    thetaEst = permThetaByAnother(params, thetaOrig, thetaEst);
    showTwoThetas(params, thetaOrig, thetaEst, true)
    showTwoThetas(params, thetaOrig, thetaEst, false)
end


function showTwoThetas(params, thetaOrig, thetaEst, withExponent)
    DOT_SIZE = 20
    thetaOrigMat = thetaToMat(params, thetaOrig);
    thetaEstMat = thetaToMat(params, thetaEst);
    inds = [params.m, params.n ^ params.order, params.k, 1]
    colors = ['b', 'r', 'g', 'm']
    figure;
    hold on;
    if withExponent
        thetaOrigMat = exp(thetaOrigMat);
        thetaEstMat = exp(thetaEstMat);
        minVal = floor(min([thetaEstMat(:); thetaOrigMat(:)]))
        maxVal = ceil(max([thetaEstMat(:); thetaOrigMat(:)]))
        plot([minVal, maxVal], [minVal, maxVal])
    end
    for i=1:4
        origMat = thetaOrigMat(:, 1:inds(i));
        thetaOrigMat = thetaOrigMat(:, inds(i)+1:end);
        estMat = thetaEstMat(:, 1:inds(i));
        thetaEstMat = thetaEstMat(:,  inds(i)+1:end);
        scatter(exp(origMat(:)), exp(estMat(:)), DOT_SIZE, colors(i), 'filled');
    end
    legend('x=y', 'T', 'E', 'G', 'startT')
    title('Learned Parameters vs True Parameters')
    xlabel('True')
    ylabel('Estimated')
end

function theta = permThetaByAnother(params, thetaOrig, thetaEst)
    perm = findCorrectThetaPermute(params, thetaOrig, thetaEst);
    theta = permTheta(thetaEst, perm);
end

function theta = permTheta(theta, perm)
    theta.T = theta.T(perm, :);
    theta.startT = theta.startT(perm);
    theta.G = theta.G(perm, :);
    for i = 1:length(perm)
        theta.E(i, :) = theta.E(perm(i), :);
    end
end

function mat = thetaToMat(params, theta)
    mat = zeros(params.m, 1 + params.m + (params.n ^ params.order) + params.k);
    for i = 1:params.m
        mat(i, :) = [theta.T(i, :), theta.E(i, :), theta.G(i, :), theta.startT(i)];
    end
end

% perm - m x 1
function perm = findCorrectThetaPermute(params, thetaOrig, thetaEst)
    % to vec
    vectorizedOrig = thetaToMat(params, thetaOrig);
    vectorizedEst = thetaToMat(params, thetaEst);
    distMat = vectorizedOrig * vectorizedEst';
    % distMat = squareform(pdist(vectorized));
    perm = zeros(params.m, 1);
    for i = 1:params.m
        [~, I] = max(distMat(:));
        [I_row, I_col] = ind2sub(size(distMat), I);
        perm(I_col) = I_row;
        distMat(I_row, :) = -1;
        distMat(:, I_col) = -1;
    end
end


function [test, train] = crossValidationSplit(params, mergedPeaksMin, testTrainRatio)
    L = size(mergedPeaksMin.seqs, 2);
    X = mergedPeaksMin.seqs;
    % X = cat(2, X, fliplr(5-X));
    % N x k x L

    pcPWMp = misc.preComputePWMp(X, params);
    N = size(X, 1);
    trainMask = true(N, 1);
    trainMask(randperm(N, floor(N * testTrainRatio))) = false;
    train.title = 'Train';
    test.title = 'Test';
    train.X = X(trainMask, :);
    test.X = X(~trainMask, :);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
    if isfield(mergedPeaksMin, 'Y')
        train.Y = mergedPeaksMin.Y(trainMask, :);
        test.Y = mergedPeaksMin.Y(~trainMask, :);
    end
    if isfield(mergedPeaksMin, 'Y2')
        train.Y2 = mergedPeaksMin.Y2(trainMask, :);
        test.Y2 = mergedPeaksMin.Y2(~trainMask, :);
    end
end

% perm - m x 1
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


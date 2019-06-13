
function mainScatterEM()
    conf.doESharing = false;
    conf.startWithBackground = false;
    conf.doBound = true;
    conf.maxIters = 1000;
    conf.canCrossLayer = true;
    conf.patience = 3;
    conf.Xpercents = [0.1, 0.5, 1];
    conf.L = 800;
    conf.N = 1000;
    conf.withExponent = false;
    conf.repeat = 1;

    conf.order = 1;
    conf.m = 3;
    conf.k = 5;
    conf.backgroundAmount = 0
    main(conf);

    conf.doBound = false;
    main(conf);

    conf.order = 2;
    conf.m = 7;
    conf.k = 10;
    conf.backgroundAmount = 1
    main(conf);

    conf.order = 3;
    conf.m = 10;
    conf.k = 30;
    conf.backgroundAmount = 2
    main(conf);
end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;
    pcPWMp = misc.preComputePWMp(mergedPeaksMin.seqs, params);
    for Xpercent = conf.Xpercents
        outpath = sprintf('m%dk%dp%do%db%d.jpg', conf.m, conf.k, floor(100 * Xpercent), conf.order, conf.doBound);
        subtitle = sprintf('m=%d, k=%d, %d%% of data', conf.m, conf.k, floor(100 * Xpercent));
        dataset.title = subtitle;
        dataset.X = mergedPeaksMin.seqs(1:floor(conf.N * Xpercent), :);
        dataset.theta = mergedPeaksMin.theta;
        dataset.pcPWMp = pcPWMp(1:floor(conf.N * Xpercent), :, :);
        [thetaEst, ~] = EM.EM(dataset, params, conf.maxIters, conf.doESharing, conf.doBound, conf.patience, conf.repeat);
        show.showTwoThetas(params, thetaOrig, thetaEst, conf.withExponent, subtitle, outpath);
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
    train.X = X(trainMask, :);
    test.title = 'Test';
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


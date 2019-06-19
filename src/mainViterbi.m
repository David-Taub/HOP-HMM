
function mainScatterEM()
    conf.doESharing = false;
    conf.startWithBackground = false;
    conf.maxIters = 1000;
    conf.canCrossLayer = true;
    conf.patience = 4;
    conf.Xpercents = [0.01, 0.25, 0.5, 0.75, 1];
    conf.L = 1000;
    conf.N = 10000;
    conf.withExponent = false;
    conf.repeat = 2;
    conf.order = 2;
    conf.m = 7;
    conf.k = 10;
    conf.backgroundAmount = 1
    conf.doBound = false;
    main(conf);

end

function main(conf)
    dbstop if error
    close all;
    params = misc.genParams(conf.m, conf.k, conf.backgroundAmount, conf.L, conf.order);
    mergedPeaksMin = mainGenSequences(conf.N, conf.L, params, conf.startWithBackground);
    thetaOrig = mergedPeaksMin.theta;
    outpath = sprintf('viterbi_m%dk%do%db%dN%dL%d.jpg', conf.m, conf.k, conf.order, conf.doBound, subN, conf.L);
    subtitle = sprintf('m=%d, k=%d, %d%% of data', conf.m, conf.k);
    dataset.title = subtitle;
    dataset.X = mergedPeaksMin.seqs;
    dataset.theta = mergedPeaksMin.theta;
    dataset.Y = mergedPeaksMin.Y;
    dataset.Y2 = mergedPeaksMin.Y2;
    dataset.pcPWMp = misc.preComputePWMp(mergedPeaksMin.seqs, params);
    [thetaEst, ~] = EM.EM(dataset, params, conf.maxIters, conf.doESharing, conf.doBound, conf.patience, conf.repeat);
    [~, ~, ~, ~, ~, psi] = EM.EStep(params, thetaEst, dataset.X, dataset.pcPWMp);
    YEstViterbi = classify(thetaEst, params, dataset);
    % N x m x k x L
    estMask = genPWMMask(params, YEstViterbi, conf.N, conf.L);
    trueMask = genPWMMask(params, YEstViterbi, conf.N, conf.L);

    estPosVals = psi(estMask);
    estNegVals = psi(~estMask);

    truePosVals = psi(trueMask);
    trueNegVals = psi(~trueMask);

    show.distributionPlot({truePosVals, estPosVals}, 'histOri', 'right', 'color', 'r', 'widthDiv', [2 2], 'showMM', 0);
    show.distributionPlot({trueNegVals, estNegVals}, 'histOri', 'left', 'color', 'b', 'widthDiv', [2 1], 'showMM', 0);

    show.violin({estPosVals, estNegVals, truePosVals, trueNegVals});




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
    YEstViterbi = misc.viterbi(params, theta, dataset.X, dataset.pcPWMp);
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


% Y - N x L x 2
% mask - N x m x k x L
function mask = genPWMMask(params, Y, N, L)
    mask = [];
    for l = 1:params.k
        layerMask = [];
        for i = 1:params.m
            PWMStateMask = (Y(:, 2:end, 2) > l) && (Y(:, 1:end - 1, 2) == 0) && (Y(:, 1:end - 1, 1) == i);
            % N x L
            PWMStateMask = cat(2, PWMStateMask, false(N, 1, 1));
            layerMask = cat(3, layerMask, PWMStateMask);
        end
        assert(size(layerMask) == [N, L, params.m]);
        layerMask = permute(layerMask, [1, 3, 4, 2]);
        mask = cat(3, mask, layerMask);
    end
    assert(size(mask) == [N, params.m, params.k, L]);
end

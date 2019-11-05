%TODO: Gen Y by the beds, permute theta, show Viterbi errors over iterations
%TODO: Download h3k27acc bigwigs, print normalized acc heatmaps of regions classified by Viterbi
function mainRealData()
    conf.startWithBackground = false;
    conf.doEnhSpecific = true;
    conf.seqsPerTissue = 1000;
    conf.maxIters = 25;
    conf.repeat = 1;
    conf.canCrossLayer = true;
    conf.patience = 4;
    conf.L = 2000;
    conf.peakMaxLength = 1000;
    conf.peakMinL = 200;
    conf.peakMaxL = 1500;
    conf.withExponent = false;
    conf.order = 3;
    conf.m = 3;
    conf.k = 20;
    conf.withBackground = true;
    conf.withGenes = false;
    conf.minSamplesCount = 10;
    conf.sequencesToShow = 5;
    conf.backgroundAmount = 1;
    conf.doESharing = false;
    conf.doGTBound = true;
    conf.doResampling = false;
    conf.topPercent = 0.5;
    conf.tissueList = [2, 37];
    % conf.tissueList = [3, 23];
    % conf.tissueList = [2, 18];
    conf.startTUniform = false;
    main(conf);
end


function main(conf)
    dbstop if error
    close all;
    mergedPeaksMin = peaks.minimizeMergePeak(conf.topPercent, conf.doEnhSpecific, conf.withBackground, conf.withGenes,...
                                             conf.seqsPerTissue, conf.L, conf.peakMinL, conf.peakMaxL, conf.tissueList,...
                                             conf.minSamplesCount);
    N = size(mergedPeaksMin.seqs, 1);
    testTrainRatio = 0.01;
    selectedPWMs = misc.PWMsFeatureSelect(mergedPeaksMin, conf.k);
    params = misc.genParams(conf.m, selectedPWMs, conf.backgroundAmount, conf.L, conf.order, ...
                            conf.doESharing, conf.doGTBound, conf.doResampling);
    [test, train] = misc.crossValidationSplit(params, mergedPeaksMin, testTrainRatio);
    [theta, ~] = EM.EM(train, params, conf.maxIters, conf.patience, conf.repeat, conf.startTUniform);
    show.showTheta(theta);
    outpath = sprintf('real_posterior_m%da%dk%do%db%dN%dL%d.jpg', conf.m, conf.backgroundAmount, ...
                      conf.k, conf.order, conf.doGTBound, N, conf.L);
    show.seqSampleCertaintyReal(params, theta, train, outpath, mergedPeaksMin.tissueEIDs);
    % seqSampleCertaintyReal(params, theta, test, conf.sequencesToShow, outpath);

    % YEst = misc.viterbi(params, theta, train.X, train.pcPWMp);
    % theta = permuteTheta(theta, params, train.Y(:, :), YEst(:, :, 1));

    % classify(theta, params, train);
    % [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, train.X, train.pcPWMp);
    % show.seqSampleCertainty(params, train.Y, gamma, psi, 8, false);

    % classify(theta, params, test);
    % [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, test.X, test.pcPWMp);
    % show.seqSampleCertainty(params, test.Y, gamma, psi, 8, false);
    keyboard
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


% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% delete(fullfile('..', 'data', 'dummyDNA.mat'));
% mainGenSequences(100, 10000, 5, 25, true, true);
% download_and_process_all.sh
% peaks.beds2matsNoSeq()
% peaks.mergePeakFiles()

% mergedPeaks = load('../data/peaks/mergedPeaks.mat', 'mergedPeaks');
% mergedPeaks = mergedPeaks.mergedPeaks;

% load(fullfile('..', 'data', 'dummyDNA.mat'), 'superEnhancers');
% superEnhancers = peaks.superEnhancerCaller(mergedPeaks, 10000);
% mainRealData(superEnhancers, 5, 40);
function mainRealData(mergedPeaksMin, m, k)
    dbstop if error
    close all;
    params = misc.genParams(m, k);
    params.NperTissue = 1000;
    testTrainRatio = 0.15;
    [test, train] = preprocess(params, mergedPeaksMin, testTrainRatio);
    % show.showTheta(mergedPeaksMin.theta);
    ttt = 0.01;
    % mergedPeaksMin.theta.G = log(exp(mergedPeaksMin.theta.G) + (ttt./params.k));
    % mergedPeaksMin.theta.T = log(exp(mergedPeaksMin.theta.T) - eye(params.m).*ttt);
    [theta, ~] = EM.EM(train.X, params, train.pcPWMp, 10);
    show.showTheta(theta);
    YEst = classify(theta, params, train);

    theta = permuteTheta(theta, params, train.Y, YEst(:, :, 1))
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, train.X, train.pcPWMp);
    show.seqSampleCertainty(params, train.Y, gamma, psi);

    classify(theta, params, test);
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, test.X, test.pcPWMp);
    show.seqSampleCertainty(params, test.Y, gamma, psi);
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

% function seqClassification(params, theta, posterior, Y)
%     % % N x m
%     % overlapsEst = matUtils.logMatSum(posterior, 3);
%     % % N x m
%     % aucRocs = zeros(params.m, 1);
%     % for i = 1:params.m
%     %     aucRocs(i) = matUtils.getAucRoc(overlapsEst(Y(:, 1) == i, i),overlapsEst(Y(:, 1) ~= i, i), false, false);
%     % end
%     % aucRocs
%     % [~, YEstMax] = max(posterior, [], 2);
%     % YEstMax = permute(YEstMax, [1, 3, 2]);
%     % accuricy = sum(sum(YEstMax == Y, 1), 2) ./ (N * L)
%     % fitcnb(overlapsEst, Y)
%     % [~, YEstMax] = max(posterior - repmat(mean(posterior, 1), [N, 1]), [], 2);
%     % YEstMax = permute(YEstMax, [1, 3, 2]);
%     % accuricy = sum(sum(YEstMax == Y, 1), 2) ./ (N * L)


%     subplot(1,5,1);
%     imagesc(theta.G); colorbar;title('G')
%     subplot(1,5,2);
%     imagesc(theta.E(:,:)); colorbar;title('E')
%     subplot(1,5,3);
%     imagesc(YSeqEstMax1hot'); colorbar;
%     title('posterior of first Seq')
%     subplot(1,5,4);
%     imagesc(overlapsEst); colorbar;
%     title('Y estimation')
%     subplot(1,5,5);
%     imagesc(matUtils.vec2mat(Y', params.m)');  colorbar;
%     title('Y real')
%     %%%%%%%%%%%%%%
% end

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

% A - l x 1
function newMask = takeTopN(A, n, mask)
    A(~mask) = -inf;
    [~, inds] = sort(A, 1);
    newMask = false(size(A, 1), 1);
    newMask(inds(end-n+1:end)) = true;
end

function [test, train] = preprocess(params, mergedPeaksMin, testTrainRatio)
    L = size(mergedPeaksMin.seqs, 2);
    X = mergedPeaksMin.seqs;
    % X = cat(2, X, fliplr(5-X));
    % N x k x L
    % k x n x J

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
        mergedPeaksMin.Y2
        train.Y2 = mergedPeaksMin.Y2(trainMask, :);
        test.Y2 = mergedPeaksMin.Y2(~trainMask, :);
    end
end


function theta = catThetas(params, thetas)
    thetas = [thetas{:}];
    params.m = length(thetas);
    theta = misc.genThetaUni(params);
    theta.G = reshape([thetas.G], [params.k, params.m])';
    %theta.F = [thetas.F]';
    theta.T = eye(params.m) * (1 - (params.m * params.tEpsilon)) + params.tEpsilon;
    theta.T = log(theta.T .* repmat(1-sum(exp(theta.G), 2), [1, params.m]));
    theta.E = zeros([params.m, ones(1, params.order) * params.n]);
    for i = 1:params.m
        theta.E(i, :) = thetas(i).E(:);
    end
end

% P - u x 1
% Q - u x 1
function ret = relativeEntropy(P, Q)
    ret = sum(P .* (log(P) - log(Q)), 1);
end

% % P - u x 1
% % Q - u x 1
% function ret = mse(P, Q)
%     ret = mean((P - Q) .^ 2, 1);
% end
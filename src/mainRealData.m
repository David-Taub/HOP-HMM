% mainGenSequences();
% load(fullfile('..', 'data', 'dummyDNA.mat'));
% mainPWM(pcPWMp, X, Y);
% mergedPeaksMin = mainGenSequences(1000, 600, 2, true);
% mergedPeaksMin = load(fullfile('..', 'data', 'dummyDNA.mat'));
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% mainGenSequences(100, 10000, 5, 25, true, true);
% mergedPeaksMin = load(fullfile('..', 'data', 'dummyDNA.mat'));
% mainRealData(mergedPeaksMin, 5, 25);
% clear all; mainGenSequences(5000, 500, 3, 6, true, true);  mainRealData(mergedPeaksMin, 5);

function mainRealData(mergedPeaksMin, m, k)
    dbstop if error
    close all;
    params = misc.genParams(m, k);
    params.NperTissue = 1000;
    testTrainRatio = 0.15;
    [test, train] = preprocess(params, mergedPeaksMin, testTrainRatio);
    [theta, ~] = EM.EM(train.X, params, train.pcPWMp, 400);
    show.showTheta(theta);
    classify(theta, params, train);
    classify(theta, params, test);
    keyboard
end


% Y - N x L
function loss = classify(theta, params, dataset)
    [N, L] = size(dataset.X);
    % N x m x L
    [alpha, beta, pX, xi, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
    % EM.drawStatus(theta, params, gamma);
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    % N x m x L
    YOneHot = permute(matUtils.mat23Dmat(dataset.Y, params.m), [1, 3, 2]);
    % N x L
    certainty = reshape(posterior(YOneHot), [N, L]);
    loss = mean(log(1-certainty(:)));
    YEstViterbi = misc.viterbi(theta, params, dataset.X, dataset.pcPWMp);
    YEstViterbiAcc = YEstViterbi(:, :, 1) == dataset.Y; YEstViterbiAcc = sum(YEstViterbiAcc(:)) ./ length(YEstViterbiAcc(:));

    YEstMax = maxPostEstimator(theta, params, psi, gamma);
    YEstMaxAcc = YEstMax(:, :, 1) == dataset.Y; YEstMaxAcc = sum(YEstMaxAcc(:)) ./ length(YEstMaxAcc(:));

    fprintf('Avg log loss: %.2f\n', loss);
    show.seqSampleCertainty(params, dataset.Y, certainty);
    figure
    subplot(1,6,1);
    imagesc(YEstViterbi(:,:,1)); colorbar; title(['Viterbi States', num2str(YEstViterbiAcc)]);
    subplot(1,6,2);
    imagesc(YEstViterbi(:,:,2)); colorbar;title('Viterbi Motifs');
    subplot(1,6,3);
    imagesc(YEstMax(:,:,1)); colorbar; title(['MaxPosterior States', num2str(YEstMaxAcc)]);
    subplot(1,6,4);
    imagesc(YEstMax(:,:,2)); colorbar;title('MaxPosterior Motifs');
    subplot(1,6,5);
    subplot(1,6,5);
    imagesc(dataset.Y); colorbar;title(['Real States (', dataset.title, ')']);
    if isfield(dataset, 'Y2')
        subplot(1,6,6);
        imagesc(dataset.Y2); colorbar;title('Real Motifs');
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
    trainMask = rand(N, 1) > testTrainRatio;
    train.title = 'Train';
    test.title = 'Test';
    train.X = X(trainMask, :);
    test.X = X(~trainMask, :);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
    if isfield(mergedPeaksMin, 'Y2')
        Y2 = mergedPeaksMin.Y2;
        train.Y = Y2(trainMask, :, 1);
        test.Y = Y2(~trainMask, :, 1);
        train.Y2 = Y2(trainMask, :, 2);
        test.Y2 = Y2(~trainMask, :, 2);
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
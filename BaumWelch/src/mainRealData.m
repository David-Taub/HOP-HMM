% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% mainPWM(pcPWMp, X, Y);
% mergedPeaksMin = mainGenSequences(1000, 600, 2, true);
% mergedPeaksMin = load(fullfile('data', 'dummyDNA.mat'));
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
% mainRealData(mergedPeaksMin);

function mainRealData(mergedPeaksMin)
    dbstop if error
    close all;
    params = misc.genParams();
    params.NperTissue = 500;
    testTrainRatio = 0.15;
    r = size(mergedPeaksMin.overlaps, 2);
    output = zeros(r, r, 2);

    % iterating over all 2 tissues pairs, each iteration around 500 sequences are picked from each tissue. then 90% of the
    % sequences are trained, then the thetas are merged and the 10% are classified using the merged thetas (with 2 floors). for
    % each classification, a rocAuc is calculated and shown in imagesc
    thetas = {};
    G = zeros(r, params.k);
    % delete(fullfile('data', 'precomputation', 'pcPWMp.mat'));
    % [test, train] = preprocess(params, mergedPeaksMin, testTrainRatio, [3, 14, 15, 9, 2]);
    [test, train] = preprocess(params, mergedPeaksMin, testTrainRatio, [1:3]);
    m = max(train.Y(:));
    for i = 1 : m
        X = train.X(train.Y(:, 1)==i, :);
        params.m = 1;
        pcPWMp = train.pcPWMp(train.Y(:, 1)==i, :, :);
        % N x k x L
        theta = learnSingleMode(X, params, pcPWMp, 7);
        thetas{i} = theta;
        drawnow;
    end
    params.m = m;
    theta = catThetas(params, thetas);
    classify(theta, params, train.X, train.pcPWMp, train.Y)
    classify(theta, params, test.X, test.pcPWMp, test.Y)
end


% Y - N x L
function loss = classify(theta, params, X, pcPWMp, Y)
    [N, L] = size(X);
    % N x m x L
    [alpha, beta, pX, xi, gamma, psi] = EM.EStep(params, theta, X, pcPWMp);
    % EM.drawStatus(theta, params, gamma);
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    % N x m x L
    YOneHot = permute(matUtils.mat23Dmat(Y, params.m), [1, 3, 2]);
    % N x L
    certainty = reshape(posterior(YOneHot), [N, L]);
    loss = mean(log(1-certainty(:)))
    keyboard
    YEst = misc.viterbi(theta, params, X, pcPWMp);
    subplot(1,3,1);
    imagesc(YEst(:,:,1)); colorbar;
    subplot(1,3,2);
    imagesc(YEst(:,:,2)); colorbar;
    subplot(1,3,3);
    imagesc(Y); colorbar;

    % subplot(1,3,1);imagesc(O1); colorbar;subplot(1,3,2);imagesc(O2(:,:,2)); colorbar;subplot(1,3,3);imagesc(O2(:,:,1)); colorbar;
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

%     %%%%%%%%%%%%%%
%     % Figures
%     %%%%%%%%%%%%%%
%     keyboard
%     sequencesToShow = 10;
%     figure;
%     inds = randsample(N, sequencesToShow);
%     for i = 1:sequencesToShow
%         subplot(sequencesToShow, 1, i);
%         YOneHot = matUtils.vec2mat(Y(inds(i), :), params.m);
%         repCertainty = YOneHot .* repmat(certainty(inds(i), :), [params.m, 1]);
%         repCertainty = [repCertainty; (1-repCertainty) .* double(repCertainty~=0)];
%         bar(repCertainty', 'stacked')
%     end
%     figure;
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

function [test, train] = preprocess(params, mergedPeaksMin, testTrainRatio, tissueIds)
    L = size(mergedPeaksMin.seqs, 2);
    overlaps = mergedPeaksMin.overlaps;
    mask = true(size(overlaps, 1), 1);
    fprintf('%d (%.2f) - > ', sum(mask, 1), sum(mask, 1)/size(overlaps, 1))
    % mask = mask & mergedPeaksMin.lengths >= 0.3*L;
    % fprintf('%d (%.2f) - > ', sum(mask, 1), sum(mask, 1)/size(overlaps, 1))
    mask = mask & mergedPeaksMin.lengths <= 3*L;
    fprintf('%d (%.2f) - > ', sum(mask, 1), sum(mask, 1)/size(overlaps, 1))
    mask = mask & (sum(overlaps > 0, 2) <= 2);
    fprintf('%d (%.2f) - > ', sum(mask, 1), sum(mask, 1)/size(overlaps, 1))
    mask = mask & (sum(overlaps(:, tissueIds) > 0, 2) == 1);
    fprintf('%d (%.2f) - > ', sum(mask, 1), sum(mask, 1)/size(overlaps, 1))
    % mask = mask & mod(1:size(mask,1), 3).' == 0;
    N = min([params.NperTissue, sum(overlaps(mask, tissueIds)>0, 1)]);
    maskTop = false(size(mask, 1), 1);
    for id = tissueIds
        maskTop = maskTop | takeTopN(overlaps(:, id), N, mask);
    end
    mask = mask & maskTop;
    fprintf('%d (%.2f)\n', sum(mask, 1), sum(mask, 1)/size(overlaps, 1))

    overlaps = overlaps(mask, :);
    overlaps = overlaps(:, tissueIds);
    X = mergedPeaksMin.seqs(mask, :);
    [overlaps, seqInd] = sortrows(overlaps);
    X = X(seqInd, :);
    Y = (overlaps > 0) * (1:size(overlaps, 2))';
    % X = cat(2, X, fliplr(5-X));
    % N x k x L
    % k x n x J

    pcPWMp = misc.preComputePWMp(X, params);
    N = size(X, 1);
    trainMask = rand(N, 1) > testTrainRatio;
    train.X = X(trainMask, :);
    train.Y = repmat(Y(trainMask), [1, L]);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.X = X(~trainMask, :);
    test.Y = repmat(Y(~trainMask), [1, L]);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
end


% pcPWMp - N x k x L
% X - N x L

function [theta] = learnSingleMode(X, params, pcPWMp, maxIter)
    [theta, ~] = EM.EMJ(X, params, pcPWMp, maxIter);
end

function theta = meanMergeTheta(params, thetas)
    thetas = [thetas{:}];
    parts = length(thetas);
    theta = misc.genThetaUni(params);
    theta.G = log(mean(reshape(exp([thetas.G]), [params.m, params.k, parts]), 3));
    theta.F = log(mean(exp([thetas.F]), 2));
    theta.T = log(mean(reshape(exp([thetas.T]), [params.m, params.m, parts]), 3));
    theta.E = zeros([params.m, ones(1, params.order) * params.n]);
    for i = 1:parts
        theta.E = theta.E + exp(thetas(i).E);
    end
    theta.E = log(theta.E ./ parts);
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
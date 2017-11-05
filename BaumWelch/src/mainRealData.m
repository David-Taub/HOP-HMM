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
        theta = learnSingleMode(X, params, pcPWMp, 3);
        thetas{i} = theta;
        drawnow;
    end
    params.m = m;
    theta = catThetas(params, thetas);
    classify(theta, params, train.X, train.pcPWMp, train.Y)
    classify(theta, params, test.X, test.pcPWMp, test.Y)
end



% gamma - N x m x L
% psi - N x m x k x L
function Yest = genEstimation(params, theta, gamma, psi)
    [N, ~, L] = size(gamma);
    % N x L x m
    gammaPer = permute(gamma, [1, 3, 2]);
    [gammaMaxVals, YestBaseStates] = max(gammaPer, [], 3);

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
    YestBaseStates(subStateMask) = baseStates(subStateMask);
    YestSubStates = zeros(N, L);
    YestSubStates(subStateMask) = subStates(subStateMask);
    Yest = cat(3, YestBaseStates, YestSubStates);

    % Yest1 = max(matUtils.logAdd(matUtils.logMatSum(skewedPsi, 4), gammaPer), [], 3);

end
% Y - N x L
function loss = classify(theta, params, X, pcPWMp, Y)
    [N, L] = size(X);
    fprintf('Calculating alpha...\n')
    % N x m x L
    alpha = EM.forwardAlgJ(X, theta, params, pcPWMp);
    fprintf('Calculating beta...\n')
    beta = EM.backwardAlgJ(X, theta, params, pcPWMp);
    % N x 1
    pX = EM.makePx(alpha, beta);
    fprintf('Calculating Gamma...\n')
    % N x m x L
    gamma = EM.makeGamma(params, alpha, beta, pX);
    % N x m x k x L
    psi = EM.makePsi(alpha, beta, X, params, theta, pcPWMp, pX);
    % EM.drawStatus(theta, params, gamma);
    % N x m x L
    YEst = calcPosterior(params, gamma, psi, N);
    % N x m x L
    YOneHot = permute(matUtils.mat23Dmat(Y, params.m), [1, 3, 2]);
    % N x L
    certainty = reshape(YEst(YOneHot), [N, L]);
    loss = mean(log(1-certainty(:)))

    % N x m
    YSeqEst = matUtils.logMatSum(YEst, 3);
    % N x 1
    [~, YSeqEstMax] = max(YSeqEst, [], 2);
    % N x m
    aucRocs = zeros(params.m, 1);
    for i = 1:params.m
        aucRocs(i) = matUtils.getAucRoc(YSeqEst(Y(:, 1) == i, i),YSeqEst(Y(:, 1) ~= i, i), false, false);
    end
    aucRocs
    [~, YEstMax] = max(YEst, [], 2);
    YEstMax = permute(YEstMax, [1, 3, 2]);

    accuracy = sum(sum(YEstMax == Y, 1), 2) ./ (N * L)

    YSeqEstMax1hot = matUtils.vec2mat(YSeqEstMax, params.m);

    %%%%%%%%%%%%%%
    % Figures
    %%%%%%%%%%%%%%
    keyboard
    sequencesToShow = 10;
    figure;
    inds = randsample(N, sequencesToShow);
    errors = repmat([1, 0.3, 0.3], [params.m, 1]);
    colormap([winter(params.m);errors])
    for i = 1:sequencesToShow
        subplot(sequencesToShow, 1, i);
        YOneHot = matUtils.vec2mat(Y(inds(i), :), params.m);
        repCertainty = YOneHot .* repmat(certainty(inds(i), :), [params.m, 1]);
        repCertainty = [repCertainty; (1-repCertainty) .* double(repCertainty~=0)];
        bar(repCertainty', 'stacked')
        xlim([1, L])
        if i < sequencesToShow
            set(gca,'xtick', [])
        end
        if i == 1
            title(sprintf('Posterior of Correct Floor\nLetter accuracy: %.2f\nLog loss: %.2f', accuracy, loss))
        end
        set(gca,'ytick', [])
        ylabel(['Seq ', num2str(inds(i))]);
    end
    xlabel('Position');
    l = strsplit(num2str(1:params.m));
    l = strcat({'Tissue type '}, l);
    l{params.m+1} = 'error';
    legend(l);


    figure;
    subplot(1,5,1);
    imagesc(theta.G); colorbar;title('G')
    subplot(1,5,2);
    imagesc(theta.E(:,:)); colorbar;title('E')
    subplot(1,5,3);
    imagesc(YSeqEstMax1hot'); colorbar;
    title('posterior of first Seq')
    subplot(1,5,4);
    imagesc(YSeqEst); colorbar;
    title('Y estimation')
    subplot(1,5,5);
    imagesc(matUtils.vec2mat(Y', params.m)');  colorbar;
    title('Y real')
    %%%%%%%%%%%%%%

end
% posterior - N x m x L
function posterior = calcPosterior(params, gamma, psi, N)
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
    theta.T = log(eye(params.m) * (1 - (params.m * params.tEpsilon)) + params.tEpsilon);
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
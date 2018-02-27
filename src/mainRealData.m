% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% delete(fullfile('..', 'data', 'dummyDNA.mat'));
% mainGenSequences(100, 10000, 5, 25, true, true);
% download_and_process_all.sh
% peaks.beds2matsNoSeq()
% peaks.mergePeakFiles()

% mergedPeaks = load('../data/peaks/mergedPeaks.mat', 'mergedPeaks');
% mergedPeaks = mergedPeaks.mergedPeaks;

% load(fullfile('..', 'data', 'dummyDNA.mat'), 'superEnhancers');
% superEnhancers = peaks.superEnhancerCaller(mergedPeaks, 10000, 10,);
% mainRealData(superEnhancers, 5, 40);
function mainRealData(mergedPeaksMin, m, k)
    dbstop if error
    close all;
    params = misc.genParams(m, k);
    params.NperTissue = 1000;
    testTrainRatio = 0.15;
    [test, train] = preprocess(params, mergedPeaksMin, testTrainRatio);
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
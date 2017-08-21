% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% pcPWMp = BaumWelchPWM.preComputePWMp(X);
% mainPWM(pcPWMp, X, Y);
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
% mainPWM(mergedPeaksMin);


% pcPWMp - N x k x L-1
function mainPWM(mergedPeaksMin)
    [Ys, Xs, pcPWMp] = genData(mergedPeaksMin);
    params.m = 5;
    params.order = 3;
    params.n = max(Xs(:));
    [params.N, params.L] = size(Xs);
    params.J = size(BaumWelchPWM.PWMs(), 3);
    params.k = size(pcPWMp, 2);
    params.tEpsilon = 0.0005;

    [theta, ~] = learn(Xs, params, pcPWMp);
    % N x m x L + J
    [~, YsEst] = max(theta.gamma(:,:,1:end-params.J), [], 10);
    YsEst = permute(YsEst, [1,3,2]);
    calcError(Ys(:)', YsEst(:)');
end

function [Ys, Xs, pcPWMp] = genData(mergedPeaksMin)
    L = size(mergedPeaksMin.seqs, 2);
    % overlaps = mergedPeaksMin.overlaps(:, :);
    overlaps = mergedPeaksMin.overlaps(:, [1, 2]);
    mask = mergedPeaksMin.lengths >= L;
    mask = mask & (sum(overlaps > 0, 2) == 1);
    mask = mask & mod(1:size(mask,1), 15).' == 0;
    overlaps = overlaps(mask, :);
    Xs = mergedPeaksMin.seqs(mask, :);
    Ys = (overlaps(:, 1) > 0) + 1;
    [overlaps, seqInd] = sortrows(overlaps);
    Xs = Xs(seqInd, :);

    % Xs = cat(2, Xs, fliplr(5-Xs));
    fprintf('Calculating PWMs LogLikelihood\n')
    size(Xs)
    pcPWMp = BaumWelchPWM.preComputePWMp(Xs);
end


% pcPWMp - N x k x L-1+J
% XTrain - N x L
function [theta, likelihood] = learn(Xs, params, pcPWMp)
    maxIter = 20;
    [theta, likelihood] = BaumWelchPWM.EMJ(Xs, params, pcPWMp, maxIter);
end


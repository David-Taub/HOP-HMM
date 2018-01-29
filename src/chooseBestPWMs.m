% mainGenSequences();
% pcPWMp = misc.preComputePWMp(X);
% mainPWM(pcPWMp, X, Y);
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% load('../data/peaks/mergedPeaksMinimized.mat');
% mergedPeaksMin = mainGenSequences(500, 600, 2, true);
% chooseBestPWMs(mergedPeaksMin, [10, 3]);


function chooseBestPWMs(mergedPeaksMin, tissueList)
    dbstop if error
    close all;
    k = 519;
    m = length(tissueList);
    params = misc.genParams(m, k);
    % params.tEpsilon = 0;
    params.batchSize = 2;

    [dataset] = preprocess(params, mergedPeaksMin, tissueList);
    % TODO: add background to the sequences
    aucRocs = zeros(params.m, params.k);
    for tissueID = 1:params.m
        tissueMask = dataset.Y(:,1) == tissueID;
        for pwmID = 1:params.k
            % N x k x L
            pos = dataset.pcPWMp(tissueMask, pwmID, :);
            neg = dataset.pcPWMp(~tissueMask, pwmID, :);
            aucRocs(tissueID, pwmID) = matUtils.getAucRoc(pos(:), neg(:), false, true);
        end
    end
    [aucRocsSorted, PWMRank] = sort(aucRocs, 2, 'descend');
    keyboard
end



function [dataset] = preprocess(params, mergedPeaksMin, tissueList)
    L = size(mergedPeaksMin.seqs, 2);
    overlaps = mergedPeaksMin.overlaps(:, :);
    % overlaps = mergedPeaksMin.overlaps(:, tissueList);
    mask = mergedPeaksMin.lengths >= L;
    mask = mask & (sum(overlaps > 0, 2) == 1);
    % mask = mask & mergedPeaksMin.Y(:,1,1) == 1;
    mask = mask & (sum(overlaps(:, tissueList) > 0, 2) == 1);
    % mask = mask & mod(1:size(mask,1), 15).' == 0;
    overlaps = overlaps(mask, :);
    overlaps = overlaps(:, tissueList);
    X = mergedPeaksMin.seqs(mask, :);
    % Y = mergedPeaksMin.Y(mask, :);
    Y = repmat((overlaps > 0) * (1:size(overlaps, 2))', [1, L]);
    [overlaps, seqInd] = sortrows(overlaps);
    X = X(seqInd, :);
    Y = Y(seqInd, :);

    % X = cat(2, X, fliplr(5-X));
    % N x k x L
    pcPWMp = misc.preComputePWMp(X, params);
    N = size(X, 1);
    dataset.title = 'dataset';
    % dataset.X = X;
    dataset.pcPWMp = pcPWMp;
    dataset.Y = Y;
    if isfield(mergedPeaksMin, 'Y')
    end
    if isfield(mergedPeaksMin, 'Y2')
        mergedPeaksMin.Y2
        dataset.Y2 = mergedPeaksMin.Y2;
    end
end

function theta = catThetas(params, thetas)
    thetas = [thetas{:}];
    params.m = length(thetas);
    theta = misc.genThetaUni(params);
    theta.G = reshape([thetas.G], [params.k, params.m])';
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
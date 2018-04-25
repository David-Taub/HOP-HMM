
function selectedPWMs = chooseBestPWMs(mergedPeaksMin, tissueList)
    dbstop if error
    profile on
    close all;
    k = 519;
    BEST_PWM_TO_CHOOSE_PER_TISSUE = 25;
    m = length(tissueList);
    params = misc.genParams(m, k);
    params.L = size(mergedPeaksMin.seqs, 2);
    % params.tEpsilon = 0;
    params.batchSize = 2;

    for tissueID = tissueList
        mergedPeaksMin.tissueNames{tissueID}
    end

    dataset = preprocess(params, mergedPeaksMin, tissueList);
    % TODO: add background to the sequences

    % m x k
    aucRocs = oneVsAllAucRoc(params, dataset);
    [aucRocsSorted, aucRocsSortedInd] = sort(aucRocs, 2, 'descend');
    % m x BEST_PWM_TO_CHOOSE_PER_TISSUE
    selectedPWMs = aucRocsSortedInd(:, 1:BEST_PWM_TO_CHOOSE_PER_TISSUE);
    selectedPWMs = unique(selectedPWMs(:));
    bestPWMsAucRocs = aucRocsSorted(:, 1:BEST_PWM_TO_CHOOSE_PER_TISSUE);
    min(bestPWMsAucRocs(:))
    max(bestPWMsAucRocs(:))
    plot(sort(bestPWMsAucRocs(:)))
    save('..\data\precomputation\SelectedPWMs.mat', 'selectedPWMs', 'aucRocsSorted', 'aucRocsSortedInd', 'tissueList');
end

function aucRocs = oneVsAllAucRoc(params, dataset)
    expected_num_of_peaks_in_seq = 5
    aucRocs = zeros(params.m, params.k);
    Xs1H = matUtils.mat23Dmat(dataset.X, params.n);
    for pwmId = 1:params.k
        PWMLogLike = misc.PWMLogLikelihood(params, Xs1H, pwmId);
        for tissueID = 1:params.m
            tissueMask = dataset.Y(:,1) == tissueID;
            % N x L
            pos = misc.maxN(PWMLogLike(tissueMask, :), 2, expected_num_of_peaks_in_seq);
            neg = misc.maxN(PWMLogLike(~tissueMask, :), 2, expected_num_of_peaks_in_seq);
            % pos = downsample(pos, 20);
            % neg = downsample(pos, 20);
            aucRocs(tissueID, pwmId) = matUtils.getAucRoc(pos(:), neg(:), false, true);
            [v, i] = max(aucRocs(:));
            fprintf('Best 1vsAll AucRocs for tissue %d of PWM %d is %.2f / %.2f (%d)\n', tissueID, pwmId, aucRocs(tissueID, pwmId), v, i);
        end
    end
end


function dataset = preprocess(params, mergedPeaksMin, tissueList)
    params.L = size(mergedPeaksMin.seqs, 2);
    overlaps = mergedPeaksMin.overlaps(:, :);
    % overlaps = mergedPeaksMin.overlaps(:, tissueList);
    mask = mergedPeaksMin.lengths >= params.L;
    mask = mask & (sum(overlaps > 0, 2) == 1);
    % mask = mask & mergedPeaksMin.Y(:,1,1) == 1;
    mask = mask & (sum(overlaps(:, tissueList) > 0, 2) == 1);
    % mask = mask & mod(1:size(mask,1), 15).' == 0;
    overlaps = overlaps(mask, :);
    overlaps = overlaps(:, tissueList);
    sum(overlaps>0, 1)
    assert(all(sum(overlaps>0, 1)>0))
    X = mergedPeaksMin.seqs(mask, :);
    % Y = mergedPeaksMin.Y(mask, :);
    Y = repmat((overlaps > 0) * (1:size(overlaps, 2))', [1, params.L]);
    [overlaps, seqInd] = sortrows(overlaps);
    X = X(seqInd, :);
    Y = Y(seqInd, :);

    % X = cat(2, X, fliplr(5-X));
    % N x k x L
    % mf = misc.preComputePWMp(X, params);
    N = size(X, 1);
    dataset.title = 'dataset';
    dataset.X = X;
    % dataset.mf = mf;
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
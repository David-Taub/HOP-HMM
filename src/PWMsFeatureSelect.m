
function selectedPWMs = PWMsFeatureSelect(mergedPeaksMin, tissueList, PWMs, lengths)
    dbstop if error
    m = length(tissueList);
    delete('../data/precomputation/SelectedPWMs.mat');
    % profile on
    % close all;
    BEST_PWM_TO_CHOOSE_PER_TISSUE = 25;
    m = length(tissueList);
    L = size(mergedPeaksMin.seqs, 2);
    k = size(PWMs, 1);
    for tissueID = tissueList
        mergedPeaksMin.tissueNames{tissueID}
    end

    dataset = preprocess(mergedPeaksMin, tissueList);
    % TODO: add background to the sequences

    % m x k
    aucRocs = oneVsAllAucRoc(dataset.X, dataset.Y, m, k, PWMs, lengths);
    [aucRocsSorted, aucRocsSortedInd] = sort(aucRocs, 2, 'descend');
    % m x BEST_PWM_TO_CHOOSE_PER_TISSUE
    selectedPWMs = aucRocsSortedInd(:, 1:BEST_PWM_TO_CHOOSE_PER_TISSUE);
    selectedPWMs = unique(selectedPWMs(:));
    bestPWMsAucRocs = aucRocsSorted(:, 1:BEST_PWM_TO_CHOOSE_PER_TISSUE);
    min(bestPWMsAucRocs(:))
    max(bestPWMsAucRocs(:))
    plot(sort(bestPWMsAucRocs(:)))
    if not(isdir('../data/precomputation'))
        mkdir('../data/precomputation')
    end
    selectedPWMsFilepath = '../data/precomputation/SelectedPWMs.mat';
    
    save(selectedPWMsFilepath, 'selectedPWMs', 'aucRocsSorted', 'aucRocsSortedInd', 'tissueList');
    fprintf('Saved feature selected PWMs in %s\n', selectedPWMsFilepath);
end

function aucRocs = oneVsAllAucRoc(X, Y, m, k, PWMs, lengths)
    expected_num_of_peaks_in_seq = 2;
    n = 4;
    aucRocs = zeros(m, k);
    Xs1H = matUtils.mat23Dmat(X, n);
    for pwmId = 1:k
        PWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, Xs1H, pwmId);
        for tissueID = 1:m
            tissueMask = Y(:,1) == tissueID;
            % N x L
            pos = misc.maxN(PWMLogLike(tissueMask, :), 2, expected_num_of_peaks_in_seq);
            neg = misc.maxN(PWMLogLike(~tissueMask, :), 2, expected_num_of_peaks_in_seq);
            aucRocs(tissueID, pwmId) = matUtils.getAucRoc(pos(:), neg(:), false, true);
            [v, i] = max(aucRocs(:));
            fprintf('Best 1vsAll AucRocs for tissue %d of PWM %d is %.2f / %.2f (%d)\n', tissueID, pwmId, aucRocs(tissueID, pwmId), v, i);
        end
    end
end


function dataset = preprocess(mergedPeaksMin, tissueList)
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
    sum(overlaps>0, 1);
    assert(all(sum(overlaps>0, 1)>0))
    X = mergedPeaksMin.seqs(mask, :);
    % Y = mergedPeaksMin.Y(mask, :);
    Y = repmat((overlaps > 0) * (1:size(overlaps, 2))', [1, L]);
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
    if isfield(mergedPeaksMin, 'Y2')
        dataset.Y2 = mergedPeaksMin.Y2;
    end
end

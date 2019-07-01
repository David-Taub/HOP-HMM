% for each couple of tissue, we take sequences that are unique to the tissues and background sequences
% then we choose the PWM with the max AucRoc and keep the
function [tissueList, pwmsList] = mainSelectTissues(mergedPeaksMin, backgroundIndex, wantedM, wantedK)
    dbstop if error
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    title('tissue aucroc maximal PWM')
    numberOfTissues = size(mergedPeaksMin.overlaps, 2);
    [PWMs, lengths, names] = misc.PWMs();
    numberOfPWMs = length(lengths);
    aucs = zeros(numberOfTissues, numberOfTissues, numberOfPWMs);
    indexMat = zeros(numberOfTissues);
    for i = 1:numberOfTissues
        for j = i + 1 : numberOfTissues
            auc1 = compareTwoTypes(mergedPeaksMin, i, [j, backgroundIndex], PWMs, lengths);
            auc2 = compareTwoTypes(mergedPeaksMin, j, [i, backgroundIndex], PWMs, lengths);
            aucs(i, j, :) = max(auc1, auc2);
        end
    end
    % take max vals tissues
    tissueList = [];
    while length(tissueList) < m - 1
        [v, inds] = max(diffMat, [], 1);
        [v2, i2] = max(v, [], 2);
        tissueList = [tissueList, i2, inds(i2)]
        diffMat(inds(i2), i2) = -inf;
    end
    tissueList = [tissueList(1:m-1), backgroundIndex];
end

% ret - k x 1
function ret = compareTwoTypes(mergedPeaksMin, ind1, inds, PWMs, lengths)
    numberOfPWMs = length(lengths);
    EXPECTED_NUM_OF_PEAKS_IN_SEQ = 2;
    SEQS_PER_TISSUE = 300;
    n = 4;
    mask1 = mergedPeaksMin.overlaps(:, ind1) > 0;
    % SEQS_PER_TISSUE x L
    seqs1 = mergedPeaksMin.seqs(mask1, :);
    assert(size(seqs1, 1) > SEQS_PER_TISSUE);
    seqs1 = seqs1(1:SEQS_PER_TISSUE, :);
    seqs2 = [];
    for i = inds
        mask2 = any(mask, mergedPeaksMin.overlaps(:, i) > 0);
        seqs = mergedPeaksMin.seqs(mask2, :);
        seqs2 = [seqs2; seqs(1:SEQS_PER_TISSUE, :)];
    end

    for pwmId = 1:k
        % SEQS_PER_TISSUE x L x n
        pos1H = matUtils.mat23Dmat(seqs1, n);
        % 2*SEQS_PER_TISSUE x L x n
        neg1H = matUtils.mat23Dmat(seqs2, n);
        % SEQS_PER_TISSUE x L
        posPWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, pos1H, pwmId);
        % 2*SEQS_PER_TISSUE x L
        negPWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, neg1H, pwmId);
        % SEQS_PER_TISSUE x EXPECTED_NUM_OF_PEAKS_IN_SEQ
        pos = misc.maxN(posPWMLogLike(tissueMask, :), 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ);
        % 2*SEQS_PER_TISSUE x EXPECTED_NUM_OF_PEAKS_IN_SEQ
        neg = misc.maxN(negPWMLogLike(~tissueMask, :), 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ);
        % N x L
        [ret(pwmId), ~, ~] = misc.getAucRoc(pos(:), neg(:), false, true);
    end
end
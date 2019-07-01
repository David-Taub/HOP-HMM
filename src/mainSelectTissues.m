% for each couple of tissue, we take sequences that are unique to the tissues and background sequences
% then we choose the PWM with the max AucRoc and keep the
function tissueList = mainSelectTissues(mergedPeaksMin, backgroundIndex, m)
    dbstop if error
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    title('tissue aucroc maximal PWM')
    MAX_SEQS_PER_TISSUE = 300;
    numberOfTissues = size(mergedPeaksMin.overlaps, 2);
    diffMat = zeros(numberOfTissues);
    indexMat = zeros(numberOfTissues);
    [PWMs, lengths, names] = misc.PWMs();
    k = size(PWMs, 1);
    for i = 1:numberOfTissues
        for j = i + 1 : numberOfTissues
            auc = compareTwoTypes(mergedPeaksMin, i, [j, backgroundIndex] ...
                                  MAX_SEQS_PER_TISSUE, k, PWMs, lengths);
            [diffMat(i, j), indexMat(i, j)] = auc;
            imagesc(diffMat);
            drawnow;
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


function ret = compareTwoTypes(mergedPeaksMin, ind1, inds, maxSeqsPerTissue, k, PWMs, lengths)
    mask1 = mergedPeaksMin.overlaps(:, ind1) > 0;
    mask2 =  any(mergedPeaksMin.overlaps(:, inds) > 0, 2);
    seqs1 = mergedPeaksMin.seqs(mask1, :);
    seqs2 = mergedPeaksMin.seqs(mask2, :);
    seqsPerTissue = min([maxSeqsPerTissue, size(seqs1, 1), size(seqs2, 1)]);
    X = [seqs2(1:seqsPerTissue, :); seqs1(1:seqsPerTissue, :)];
    Y = [ones(seqsPerTissue, 1); ones(seqsPerTissue, 1) * 2];
    aucRocs = misc.oneVsAllAucRoc(X, Y, Z, k, PWMs, lengths)
    ret = max(aucRocs(:));
end
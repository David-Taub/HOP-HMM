% for each couple of tissue, we take sequences that are unique to the tissues and background sequences
then we choose the PWM with the max AucRoc and keep the
function tissueList = mainSelectTissues(mergedPeaksMin, backgroundIndex, m)
    dbstop if error
    figure;
    title('tissue aucroc maximal PWM')
    tissuesCount = size(mergedPeaksMin.overlaps, 2);
    diffMat = zeros(tissuesCount, tissuesCount);
    indexMat = zeros(tissuesCount, tissuesCount);
    [PWMs, lengths, names] = misc.PWMs();
    k = size(PWMs, 1);
    backMask = mergedPeaksMin.overlaps(:, backgroundIndex) > 0;
    sum(backMask, 1)
    for i = 1:tissuesCount
        tissueMask = mergedPeaksMin.overlaps(:, i) > 0;
        sum(tissueMask, 1)
        backSeqs = mergedPeaksMin.seqs(backMask, :);
        tissueSeqs = mergedPeaksMin.seqs(tissueMask, :);
        maxSeqsPerTissue = 300;
        maxSeqsPerTissue = min([maxSeqsPerTissue, size(backSeqs, 1)]);
        maxSeqsPerTissue = min([maxSeqsPerTissue, size(tissueSeqs, 1)]);
        X = [backSeqs(1:maxSeqsPerTissue, :);
             tissueSeqs(1:maxSeqsPerTissue, :)];
        Y = [ones(maxSeqsPerTissue, 1); ones(maxSeqsPerTissue, 1) * 2];
        m = 2;
        aucRocs = misc.oneVsAllAucRoc(X, Y, m, k, PWMs, lengths)
        [diffMat(i, i), indexMat(i, i)] = max(aucRocs(:));

        for j = i+1 : tissuesCount
            tissueMask1 = mergedPeaksMin.overlaps(:, i) > 0;
            tissueMask2 = mergedPeaksMin.overlaps(:, j) > 0;
            backSeqs = mergedPeaksMin.seqs(backMask, :);
            tissue1Seqs = mergedPeaksMin.seqs(tissueMask1, :);
            tissue2Seqs = mergedPeaksMin.seqs(tissueMask2, :);
            maxSeqsPerTissue = 300;
            maxSeqsPerTissue = min([maxSeqsPerTissue, size(backSeqs, 1)]);
            maxSeqsPerTissue = min([maxSeqsPerTissue, size(tissue1Seqs, 1)]);
            maxSeqsPerTissue = min([maxSeqsPerTissue, size(tissue2Seqs, 1)]);
            X = [backSeqs(1:maxSeqsPerTissue, :);
                 tissue1Seqs(1:maxSeqsPerTissue, :);
                 tissue2Seqs(1:maxSeqsPerTissue, :)];
            Y = [ones(maxSeqsPerTissue, 1);
                 ones(maxSeqsPerTissue, 1) * 2;
                 ones(maxSeqsPerTissue, 1) * 3];
            m = 3;
            aucRocs = misc.oneVsAllAucRoc(X, Y, m, k, PWMs, lengths)
            [diffMat(i, j), indexMat(i, j)] = max(aucRocs(:));
            imagesc(diffMat);
            drawnow;
        end
    end
    tissueList = [];
    while length(tissueList) < m-1
        [v, i] = max(diffMat, [], 1);
        [v2, i2] = max(v, [], 2);
        tissueList = [tissueList, i2, i(i2)]
        diffMat(i(i2), i2) = -inf;
    end
    tissueList = [tissueList(1:m-1), backgroundIndex];
end
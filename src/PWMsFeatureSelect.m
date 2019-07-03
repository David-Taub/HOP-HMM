
function selectedPWMs = PWMsFeatureSelect(mergedPeaksMin, wantedK)
    dbstop if error
    outFilepath = sprintf('../data/precomputation/SelectedPWMs_L%dwk%dm%dtissues%s.mat', L,...
                          wantedK, m, sprintf('%d', mergedPeaksMin.tissueList));
    if isfile(outFilepath)
        load(outFilepath);
        return;
    end
    [PWMs, lengths, ~] = misc.PWMs();
    m = size(mergedPeaksMin.overlaps, 2);
    L = size(mergedPeaksMin.seqs, 2);
    k = size(PWMs, 1);

    [X, Y] = preprocess(mergedPeaksMin);

    % m x k
    aucRocs = misc.oneVsAllAucRoc(X, Y, PWMs, lengths);
    [aucRocsSorted, aucRocsSortedInd] = sort(aucRocs, 2, 'descend');

    selectedPWMs = [];
    i = 1;
    while length(selectedPWMs) < wantedK
        aucRocsSorted(1:i)
        selectedPWMs = unique(aucRocsSortedInd(1:i));
    end
    save(outFilepath, 'selectedPWMs');
    fprintf('Saved feature selected PWMs in %s\n', outFilepath);
end



function [X, Y] = preprocess(mergedPeaksMin)
    X = mergedPeaksMin.seqs(mask, :);
    Y = (mergedPeaksMin.overlaps > 0) * (1:size(mergedPeaksMin.overlaps, 2))'
    [overlaps, seqInd] = sortrows(mergedPeaksMin.overlaps);
    X = X(seqInd, :);
    Y = Y(seqInd);
end

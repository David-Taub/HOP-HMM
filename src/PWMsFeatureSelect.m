
function selectedPWMs = PWMsFeatureSelect(mergedPeaksMin, tissueList, PWMs, lengths, wantedK)
    dbstop if error
    m = length(tissueList);
    outFilepath = sprintf('../data/precomputation/SelectedPWMs_N%dL%dk%dwk%dm%dtissues%s.mat', N, L, k,...
                          wantedK, m, sprintf('%d', tissueList));
    if isfile(outFilepath)
        load(outFilepath);
        return;
    end
    % profile on
    % close all;
    BEST_PWM_TO_CHOOSE_PER_TISSUE = 25;
    m = length(tissueList);
    L = size(mergedPeaksMin.seqs, 2);
    k = size(PWMs, 1);
    for tissueID = tissueList
        mergedPeaksMin.tissueNames{tissueID}
    end

    [X, Y] = preprocess(mergedPeaksMin, tissueList);

    % m x k
    aucRocs = misc.oneVsAllAucRoc(X, Y, PWMs, lengths);
    [aucRocsSorted, aucRocsSortedInd] = sort(aucRocs, 2, 'descend');
    % m x BEST_PWM_TO_CHOOSE_PER_TISSUE

    selectedPWMs = [];
    i = 1;
    while length(selectedPWMs) < wantedK
        aucRocsSorted(1:i)
        selectedPWMs = unique(aucRocsSortedInd(1:i));
    end
    save(outFilepath, 'selectedPWMs');
    fprintf('Saved feature selected PWMs in %s\n', outFilepath);
end



function [X, Y] = preprocess(mergedPeaksMin, tissueList)
    mask = sum(overlaps(:, tissueList) > 0, 2) == 1;
    overlaps = overlaps(mask, :);
    overlaps = overlaps(:, tissueList);
    assert(all(sum(overlaps > 0, 1) > 0))
    X = mergedPeaksMin.seqs(mask, :);
    Y = (overlaps > 0) * (1:size(overlaps, 2))'
    [overlaps, seqInd] = sortrows(overlaps);
    X = X(seqInd, :);
    Y = Y(seqInd);
end

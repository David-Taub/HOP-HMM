
% L - super enhancer size
% mergedPeaks - enhancers data
% finds the J regions of length L with the most enhancers starts in them
function multiEnhancers = multiEnhancerCaller(mergedPeaks, L, tissueList)
    close all;
    binSize = 500;
    amountToPick = 100;
    %todo: remove non unique enhancers
    multiEnhancers = findMultiEnhancers(mergedPeaks, L, binSize, amountToPick, tissueList);
    multiEnhancers = addSequences(multiEnhancers, L);
    multiEnhancers = addY(multiEnhancers, mergedPeaks, L);
    multiEnhancersFilepath = fullfile('..', 'data', 'multiEnhancers.mat')
    save(multiEnhancersFilepath, 'multiEnhancers');
    fprintf('Saved multiEnhancers in %s\n', multiEnhancersFilepath)
end


function enhancers = getEnhancersInRegion(mergedPeaks, chr, fromSuperEnh, toSuperEnh)
    mask = strcmp({mergedPeaks.chr}, chr);
    mask = mask & [mergedPeaks.peakTo] >= fromSuperEnh;
    mask = mask & [mergedPeaks.peakFrom] <= toSuperEnh;
    enhancers = mergedPeaks(mask);
end


function multiEnhancers = addY(multiEnhancers, mergedPeaks, L)
    m = length(mergedPeaks(1).overlap);
    N = length(multiEnhancers.chr);
    multiEnhancers.Y = ones(N, L) * (m + 1);
    for i = 1:N
        chr = multiEnhancers.chr{i};
        fromSuperEnh = multiEnhancers.from(i);
        toSuperEnh = fromSuperEnh + L - 1;
        enhancers = getEnhancersInRegion(mergedPeaks, chr, fromSuperEnh, toSuperEnh);
        assert(multiEnhancers.amount(i) == length(enhancers))
        for j = 1:length(enhancers)
            % Yval = (enhancers(j).overlap > 0) * [1:m]';
            Yval = find(enhancers(j).overlap > 0, 1);
            fromEnh = max(fromSuperEnh, enhancers(j).peakFrom) - fromSuperEnh + 1;
            toEnh = min(toSuperEnh, enhancers(j).peakTo) - fromSuperEnh + 1;
            multiEnhancers.Y(i, fromEnh : toEnh) = Yval;
        end
    end
    v = sort(unique(multiEnhancers.Y(:)));
    for i = 1:length(v)
        multiEnhancers.Y(multiEnhancers.Y == v(i)) = i;
    end
end


function multiEnhancers = addSequences(multiEnhancers, L)
    dict = peaks.fasta2mem();
    N = length(multiEnhancers.chr);
    multiEnhancers.seqs = zeros(N, L);
    for i = 1:N
        chr = multiEnhancers.chr{i};
        from = multiEnhancers.from(i);
        to = multiEnhancers.from(i) + L - 1;
        multiEnhancers.seqs(i, :) = dict(chr).Data(from:to);

    end
end


function multiEnhancers = findMultiEnhancers(mergedPeaks, L, binSize, amountToPick, tissueList)
    HIGH_PEAKS_RATIO = 0.9; %amount of high peaks used, 0=none 1=all
    MAX_PEAKS_LENGTH = 1200;
    [multiEnhancers.chr{1:amountToPick}] = deal('');
    multiEnhancers.amount = zeros(1, amountToPick);
    multiEnhancers.from = zeros(1, amountToPick);

    for chrName = unique({mergedPeaks.chr})
        chrSize = getChrSize(chrName);
        fprintf('~~%s (%d)\n', chrName{:}, chrSize);
        wantedMask = getChrMask(mergedPeaks, chrName);
        fprintf('%.3f < ', sum(wantedMask,1) ./ length(wantedMask));
        wantedMask = wantedMask & getUniqueMask(mergedPeaks);
        fprintf('%.3f < ', sum(wantedMask,1) ./ length(wantedMask));
        wantedMask = wantedMask & highFilterPeaks(mergedPeaks, HIGH_PEAKS_RATIO);
        fprintf('%.3f < ', sum(wantedMask,1) ./ length(wantedMask));
        wantedMask = wantedMask & getTissueListMask(mergedPeaks, tissueList);
        fprintf('%.3f < ', sum(wantedMask,1) ./ length(wantedMask));
        wantedMask = wantedMask & getShortMask(mergedPeaks, MAX_PEAKS_LENGTH);
        fprintf('%.3f\n', sum(wantedMask,1) ./ length(wantedMask));
        wantedPeaks = mergedPeaks(wantedMask);

        unwantedMask = false(length(mergedPeaks), 1);
        unwantedMask = unwantedMask | not(getTissueListMask(mergedPeaks, tissueList));
        fprintf('%.3f < ', sum(unwantedMask,1) ./ length(unwantedMask));
        unwantedMask = unwantedMask | not(getUniqueMask(mergedPeaks));
        fprintf('%.3f < ', sum(unwantedMask,1) ./ length(unwantedMask));
        % unwantedMask = unwantedMask | not(highFilterPeaks(mergedPeaks, HIGH_PEAKS_RATIO));
        % fprintf('%.3f < ', sum(unwantedMask,1) ./ length(unwantedMask));
        % unwantedMask = unwantedMask | not(getShortMask(mergedPeaks, MAX_PEAKS_LENGTH));
        fprintf('%.3f < ', sum(unwantedMask,1) ./ length(unwantedMask));
        unwantedMask = unwantedMask & getChrMask(mergedPeaks, chrName);
        fprintf('%.3f\n', sum(unwantedMask,1) ./ length(unwantedMask));
        unwantedPeaks = mergedPeaks(unwantedMask);

        edges = 0 : chrSize ./ binSize;
        wantedPeaksHist = histcounts(floor([[wantedPeaks.peakFrom], [wantedPeaks.peakTo]] ./ binSize), edges);
        unwantedPositions = [[[unwantedPeaks.peakFrom], [unwantedPeaks.peakTo]] ./ binSize , 1 + ([[unwantedPeaks.peakFrom], [unwantedPeaks.peakTo]] ./ binSize)];
        unwantedPeaksHist = histcounts(unwantedPositions, edges);

        wantedPeaksHist(unwantedPeaksHist > 0) = -inf;
        % divide in 2 since we count both the start and the end of the enhancers
        windows = NCumSum(1 + floor(L / binSize), wantedPeaksHist) ./ 2;
        multiEnhancers = updatemultiEnhancers(multiEnhancers, windows, chrName{:}, L, binSize, mergedPeaks, tissueList);

    end
end


% mask of peaks that appear in only one tissue
function mask = getChrMask(mergedPeaks, chr)
    mask = strcmp({mergedPeaks.chr}, chr)';
end


% mask of peaks that appear in only one tissue
function mask = getShortMask(mergedPeaks, maxLength)
    mask = ([mergedPeaks.peakTo] - [mergedPeaks.peakFrom] < maxLength)';
end


% mask of peaks that appear in only one tissue
function mask = getUniqueMask(mergedPeaks)
    overlaps = vertcat(mergedPeaks.overlap);
    mask = sum(overlaps > 0, 2) == 1;
end


% mask of peaks that appear only in given tissues
function mask = getTissueListMask(mergedPeaks, tissueList)
   overlaps = vertcat(mergedPeaks.overlap);
   overlaps(:, tissueList) = 0;
   mask = sum(overlaps, 2) == 0;
end


% mask of highest P percent of the peaks (in non unique, highest peaks is measured)
function mask = highFilterPeaks(mergedPeaks, P)
   overlaps = vertcat(mergedPeaks.overlap);
   overlapsMax = max(overlaps, [], 2);
   [~, inds] = sort(overlapsMax);
   mask = false(size(overlaps, 1), 1);
   mask(inds(end-round(length(inds) * P):end)) = true;

end


function out = NCumSum(n, v)
    out = zeros(1, size(v, 2));
    for i = 1 : n
        out = out + v;
        %shift left
        v = [v(2:end), [0]];
    end
end


function multiEnhancers = updatemultiEnhancers(multiEnhancers, windows, chrName, L, binSize, mergedPeaks, tissueList)
    assert(length(multiEnhancers) > 0)
    while max(windows, [], 2) > min(multiEnhancers.amount, [], 2)
        [toReplaceAmount, toReplace] = min(multiEnhancers.amount, [], 2);
        [~, newBestFrom] = max(windows, [], 2);
        newBestFrom = newBestFrom * binSize;
        newBestTo = newBestFrom + L - 1;

        enhancers = getEnhancersInRegion(mergedPeaks, chrName, newBestFrom, newBestTo);
        newAmount = length(enhancers);

        overlaps = vertcat(enhancers.overlap);
        overlaps(:, tissueList) = 0;
        containsOnlyTissueList = sum(overlaps(:), 1) == 0;

        if(newAmount > 0 && newAmount > toReplaceAmount && containsOnlyTissueList)
            fprintf('%s:%d (%d) <- %s:%d (%d)\n', multiEnhancers.chr{toReplace}, ...
                                                multiEnhancers.from(toReplace), ...
                                                toReplaceAmount, ...
                                                chrName, ...
                                                newBestFrom, ...
                                                newAmount)
            multiEnhancers.chr{toReplace} = chrName;
            multiEnhancers.amount(toReplace) = newAmount;
            multiEnhancers.from(toReplace) = newBestFrom;
        end
        unmarkBinsFrom = max(0, newBestFrom / binSize - floor(L / binSize));
        unmarkBinsTo = min(size(windows, 2), newBestFrom / binSize + floor(L / binSize));
        windows(unmarkBinsFrom:unmarkBinsTo) = -inf;
    end
end

function chrSize = getChrSize(chrName)
    FILE_PATH = '../data/peaks/raw_bed/hg19.chrom.sizes';
    fid = fopen(FILE_PATH);
    data = textscan(fid, '%s%d', 'delimiter','\t');
    lengths = data{2};
    chrSize = lengths(strcmp(data{1}, chrName));
    fclose(fid);
end


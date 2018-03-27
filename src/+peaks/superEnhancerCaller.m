
% L - super enhancer size
% mergedPeaks - enhancers data
% finds the J regions of length L with the most enhancers starts in them

function superEnhancers = superEnhancerCaller(mergedPeaks, L, tissueList)
    close all;
    binSize = 500;
    amountToPick = 100;
    %todo: remove non unique enhancers
    superEnhancers = findSuperEnhancers(mergedPeaks, L, binSize, amountToPick, tissueList);
    superEnhancers = addSequences(superEnhancers, L);
    superEnhancers = addY(superEnhancers, mergedPeaks, L);
    save(fullfile('..', 'data', 'superEnhancers.mat'), 'superEnhancers');
end


function enhancers = getEnhancersInRegion(mergedPeaks, chr, fromSuperEnh, toSuperEnh)
    mask = strcmp({mergedPeaks.chr}, chr);
    mask = mask & [mergedPeaks.peakTo] >= fromSuperEnh;
    mask = mask & [mergedPeaks.peakFrom] <= toSuperEnh;
    enhancers = mergedPeaks(mask);
end


function superEnhancers = addY(superEnhancers, mergedPeaks, L)
    m = length(mergedPeaks(1).overlap);
    N = length(superEnhancers.chr);
    superEnhancers.Y = zeros(N, L);
    for i = 1:N
        chr = superEnhancers.chr{i};
        fromSuperEnh = superEnhancers.from(i);
        toSuperEnh = superEnhancers.from(i) + L - 1;
        enhancers = getEnhancersInRegion(mergedPeaks, chr, fromSuperEnh, toSuperEnh);
        for j = 1:length(enhancers)
            Yval = (enhancers(j).overlap > 0) * [1:m]';
            fromEnh = max(fromSuperEnh, enhancers(j).peakFrom) - fromSuperEnh + 1;
            toEnh = min(toSuperEnh, enhancers(j).peakTo) - fromSuperEnh + 1;
            superEnhancers.Y(i, fromEnh : toEnh) = Yval;
        end
    end

    v = sort(unique(superEnhancers.Y(:)));
    for i = 1:length(v)
        superEnhancers.Y(superEnhancers.Y == v(i)) = i;
    end
end


function superEnhancers = addSequences(superEnhancers, L)
    dict = peaks.fasta2mem();
    N = length(superEnhancers.chr);
    superEnhancers.seqs = zeros(N, L);
    for i = 1:N
        chr = superEnhancers.chr{i};
        from = superEnhancers.from(i);
        to = superEnhancers.from(i) + L - 1;
        superEnhancers.seqs(i, :) = dict(chr).Data(from:to);

    end
end


function superEnhancers = findSuperEnhancers(mergedPeaks, L, binSize, amountToPick, tissueList)
    HIGH_PEAKS_RATIO = 0.9; %amount of high peaks used, 0=none 1=all
    MAX_PEAKS_LENGTH = 1200;
    [superEnhancers.chr{1:amountToPick}] = deal('');
    superEnhancers.amount = zeros(1, amountToPick);
    superEnhancers.from = zeros(1, amountToPick);

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
        unwantedPeaksHist = histcounts(floor([[unwantedPeaks.peakFrom], [unwantedPeaks.peakTo]] ./ binSize), edges);
        wantedPeaksHist(unwantedPeaksHist > 0) = -inf;
        windows = NCumSum(1+floor(L / binSize), wantedPeaksHist);
        superEnhancers = updateSuperEnhancers(superEnhancers, windows, chrName, L, binSize);

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


function superEnhancers = updateSuperEnhancers(superEnhancers, windows, chrName, L, binSize)
    while max(windows, [], 2) > min(superEnhancers.amount, [], 2)
        [toReplaceAmount, toReplace] = min(superEnhancers.amount, [], 2);
        [newBest, newBestfrom] = max(windows, [], 2);
        fprintf('%s:%d (%d) <- %s:%d (%d)\n', superEnhancers.chr{toReplace}, ...
                                            superEnhancers.from(toReplace), ...
                                            superEnhancers.amount(toReplace), ...
                                            chrName{:}, ...
                                            newBestfrom * binSize, ...
                                            newBest)

        superEnhancers.chr{toReplace} = chrName{:};
        superEnhancers.amount(toReplace) = newBest;
        superEnhancers.from(toReplace) = newBestfrom * binSize;
        unmarkBinsFrom = max(0,newBestfrom - floor(L / binSize));
        unmarkBinsTo = min(size(windows, 2) , newBestfrom + floor(L / binSize));
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


% mergedPeaks = load('../data/peaks/mergedPeaks.mat', 'mergedPeaks');
% superEnhancers = peaks.superEnhancerCaller(mergedPeaks, 10000);
% L - super enhancer size
% mergedPeaks - enhancers data
% finds the J regions of length L with the most enhancers starts in them
function superEnhancers = superEnhancerCaller(mergedPeaks, L)
    close all;
    binTicks = 2000;
    amountToPick = 50;
    superEnhancers = findSuperEnhancers(mergedPeaks, L, binTicks, amountToPick);
    superEnhancers = addSequences(superEnhancers, L);
end

function superEnhancers = addSequences(superEnhancers, L)
    dict = peaks.fasta2mem();
    for i = 1 : length(superEnhancers.chr)
        chr = superEnhancers.chr{i};
        from = superEnhancers.from(i);
        to = superEnhancers.from(i) + L - 1;
        superEnhancers.seq{i} = dict(chr).Data(from:to);
    end
end

function superEnhancers = findSuperEnhancers(mergedPeaks, L, binTicks, amountToPick)
    [superEnhancers.chr{1:amountToPick}] = deal('');
    superEnhancers.amount = zeros(1, amountToPick);
    superEnhancers.from = zeros(1, amountToPick);
    for chrName = unique({mergedPeaks.chr})
        chrSize = getChrSize(chrName);
        fprintf('~~%s (%d)\n', chrName{:}, chrSize);
        chrPeaks = mergedPeaks(strcmp({mergedPeaks.chr}, chrName));
        windows = zeros(1, floor(chrSize ./ binTicks));
        edges = 0 : chrSize ./ binTicks;
        chrHist = histcounts(floor([chrPeaks.peakFrom] ./ binTicks), edges);
        for i = 1 : floor(L / binTicks)
            windows = windows + chrHist;
            chrHist = [chrHist(2:end), [0]];
        end
        while max(windows, [], 2) > min(superEnhancers.amount, [], 2)
            [toReplaceAmount, toReplace] = min(superEnhancers.amount, [], 2);
            [newBest, newBestfrom] = max(windows, [], 2);
            fprintf('%s:%d (%d) <- %s:%d (%d)\n', superEnhancers.chr{toReplace}, ...
                                                superEnhancers.from(toReplace), ...
                                                superEnhancers.amount(toReplace), ...
                                                chrName{:}, ...
                                                newBestfrom * binTicks, ...
                                                newBest)
            superEnhancers.chr{toReplace} = chrName{:};
            superEnhancers.amount(toReplace) = newBest;
            superEnhancers.from(toReplace) = newBestfrom * binTicks;
            windows(max(0,newBestfrom - floor(L / binTicks))+1:newBestfrom) = -inf;
        end
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


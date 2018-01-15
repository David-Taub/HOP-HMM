% mergedPeaks = load('../data/peaks/mergedPeaks.mat', 'mergedPeaks');
% peaks.superEnhancerCaller(mergedPeaks, 10000);
% L - super enhancer size
% mergedPeaks - enhancers data
% finds the J regions of length L with the most enhancers starts in them
function superEnhancerCaller(mergedPeaks, L)
    binTicks = 1000;
    amountToPick
    for chrName = unique({mergedPeaks.chr})
        chrSize = getChrSize(chrName);
        chrPeaks = mergedPeaks(strcmp({mergedPeaks.chr}, chrName));
        windows = zeros(1, floor(chrSize ./ binTicks));
        edges = 0 : chrSize ./ binTicks;
        chrHist = histcounts(floor([chrPeaks.peakFrom] ./ binTicks), edges);
        for i = 1 : floor(L / binTicks)
            windows = windows + chrHist;
            chrHist = [chrHist(2:end), [0]];
        end
        figure
        plot(chrHist)
        drawnow;
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


% merge all peaks.mat from /cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren
% to a single mat file. These peaks thought to be enhancers, and are basically
% P300 peaks in MM9

function mergePeakFiles()
    [totalpeaks, names] = genTotalPeaks();
    mergedPeaks = genMergePeaks(totalpeaks);
    save('-v7.3', './merged/mergedPeaks.mat', 'mergedPeaks', 'names');
end

function [totalpeaks, names] = genTotalPeaks()
    totalpeaks = [];
    peaksBasePath = 'raw';
    peakFiles = dir(fullfile(peaksBasePath, '*.peaks.mat'));
    for i = 1:length(peakFiles)
        peakFiles(i).name
        peaks = load(fullfile(peaksBasePath, peakFiles(i).name));
        S = [peaks.S{:}];
        for j=1:size(S, 2)
            S(j).overlap(i) = S(j).height;
            S(j).overlap(S(j).overlap == 1) = -inf;
        end
        totalpeaks = [totalpeaks, S];
    end
end

function mergedPeaks  = genMergePeaks(totalpeaks)

    mergedPeaks = totalpeaks;
    j = 1;
    for chrName = unique({totalpeaks.chr})
        chrMask = strcmp({totalpeaks.chr}, chrName{1});
        chrPeaks = totalpeaks(chrMask);
        [~, ind] = sort([chrPeaks.from]);
        mergedPeaks(j) = chrPeaks(ind(1));
        i = 2;
        while i <= length(ind)
            newPeak = chrPeaks(ind(i));
            oldPeak = mergedPeaks(j);
            if oldPeak.to > newPeak.from
                fprintf('.')
                % merge
                oldPeak.seq = [oldPeak.seq, newPeak.seq(end-(newPeak.to-oldPeak.to) + 1:end)];
                oldPeak.to = newPeak.to;
                oldPeak.pos = round((oldPeak.from + newPeak.to)/2) ;
                oldPeak.overlap = max(oldPeak.overlap, newPeak.overlap);
                oldPeak.height = max(oldPeak.height, newPeak.height);
                oldPeak.min = min(oldPeak.min, newPeak.min);
                mergedPeaks(j) = oldPeak;

            else
                fprintf('O')
                j = j + 1;
                mergedPeaks(j) = newPeak;
            end
            if mod(i,100)==0
                fprintf(' %s, %d / %d\n', chrName{1}, i, length(ind));
            end
            i = i + 1;
        end
    end
    mergedPeaks = mergedPeaks(1:j);
end

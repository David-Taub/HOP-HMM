% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src


% peaks.beds2mats(500)
% peaks.mergePeakFiles()
% load('../data/peaks/mergedPeaks.mat')
% peaks.minimizeMergePeak(mergedPeaks, 500, tissueNames);
% mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized.mat')
% chooseBestPWMs(mergedPeaksMin, [1, 2, 3, 45, 46]);

function minimizeMergePeak(mergedPeaks, L, tissueNames)
    TOP_PEAKS_HEIGHT_PERCENT_KEPT = 0.30;
    numerOfTissues = length(mergedPeaks(1).overlap);


    % mergedPeaks = removeNonDistal(mergedPeaks);
    [overlaps, lengths] = extractOverlaps(mergedPeaks, numerOfTissues);
    size(overlaps)
    [mergedPeaks, overlaps, lengths] = removeLowHeight(mergedPeaks, overlaps, lengths, numerOfTissues, TOP_PEAKS_HEIGHT_PERCENT_KEPT);
    size(overlaps)
    seqs = trimSeqs(mergedPeaks, L);
    [seqs, overlaps, lengths] = removeNonLetters(seqs, overlaps, lengths);
    size(overlaps)
    lengths = double(lengths);

    fprintf('save\n');
    save('../data/peaks/mergedPeaksMinimized.mat', '-v7.3', 'seqs', 'overlaps', 'lengths', 'tissueNames')
end


function [seqs, overlaps, lengths] = removeNonLetters(seqs, overlaps, lengths)
    fprintf('non letters\n');
    mask = max(seqs, [], 2) <= 4;
    seqs = seqs(mask, :);
    overlaps = overlaps(mask, :);
    lengths = lengths(mask);
end

function seqs = trimSeqs(mergedPeaks, L)

    fprintf('sequence\n');
    seqsCells = {mergedPeaks.seq};
    seqs = zeros(length(mergedPeaks), L);
    for i = 1:length(seqsCells)
        seq = seqsCells{i};
        center = round(length(seq)/2);
        % seqs(i, :) = nt2int(seq(center-L/2+1:center+L/2));
        start_pos = center - L / 2 +1;
        end_pos = start_pos + L - 1 ;
        seqs(i, :) = seq(start_pos:end_pos);
    end
end

function [mergedPeaks, overlaps, lengths] = removeLowHeight(mergedPeaks, overlaps, lengths, numerOfTissues, topPercent)
    fprintf('height\n');
    % remove low peaks
    mask = false(length(mergedPeaks), 1);
    for i = 1:numerOfTissues
        [vals, ind] = sort(overlaps(:,i), 'descend');
        amountToKeep = round(sum(overlaps(:,i) > 0)*topPercent);
        if vals(1) == vals(end)
            % not an enhancer height peak, either background or gene peak
            amountToKeep = sum(overlaps(:,i) > 0);
        end
        mask(ind(1:amountToKeep)) = true;
    end
    overlaps = overlaps(mask, :);
    mergedPeaks = mergedPeaks(mask);
    lengths = lengths(mask);
end

function [overlaps, lengths] = extractOverlaps(mergedPeaks, numerOfTissues)

    fprintf('overlaps\n');
    overlapsFlat = [mergedPeaks.overlap];
    lengths = [mergedPeaks.peakLength]';
    % N x numerOfTissues
    overlaps = reshape(overlapsFlat, [numerOfTissues, length(mergedPeaks)])';
end

function mergedPeaks = removeNonDistal(mergedPeaks)
    fprintf('Remove non distal\n');
    mergedPeaks = mergedPeaks(strcmp({mergedPeaks.class}, 'Distal'));
end
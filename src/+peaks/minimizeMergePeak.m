% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src

function minimizeMergePeak(mergedPeaks, L, tissueNames)
    TOP_PEAKS_HEIGHT_PERCENT_KEPT = 0.50;
    numerOfTissues = length(mergedPeaks(1).overlap);


    % mergedPeaks = removeNonDistal(mergedPeaks);
    [overlaps, lengths] = extractOverlaps(mergedPeaks, numerOfTissues);
    size(overlaps, 1)
    seqs = trimSeqs(mergedPeaks, L);
    size(overlaps, 1)
    [seqs, overlaps, lengths] = removeLowHeight(seqs, overlaps, lengths, numerOfTissues, TOP_PEAKS_HEIGHT_PERCENT_KEPT);
    size(overlaps, 1)
    [seqs, overlaps, lengths] = removeNonLetters(seqs, overlaps, lengths);
    size(overlaps, 1)
    [seqs, overlaps, lengths] = balanceOverlaps(seqs, overlaps, lengths);
    outFilepath = '../data/peaks/mergedPeaksMinimized.mat';
    save(outFilepath, '-v7.3', 'seqs', 'overlaps', 'lengths', 'tissueNames');
    fprintf('Saved peaks in %s\n', outFilepath);
end

function [seqs, overlaps, lengths] = balanceOverlaps(seqs, overlaps, lengths)
    MAX_SKEW_FACTOR = 3;
    counts = sum(overlaps(sum(overlaps > 0, 2) == 1, :) > 0, 1)
    shouldDecrease = counts > (median(counts, 2) * MAX_SKEW_FACTOR);
    factors = (median(counts, 2) * MAX_SKEW_FACTOR) ./ counts;
    for i = find(shouldDecrease)
        mask = sum(overlaps>0, 2) == 1;
        mask = mask & overlaps(:, i) > 0;
        mask = mask & rand(size(mask)) > factors(i);
        mask = ~mask;
        [seqs, overlaps, lengths] = reduceData(mask, seqs, overlaps, lengths);
    end
end

function [seqs, overlaps, lengths] = reduceData(mask, seqs, overlaps, lengths)
    seqs = seqs(mask, :);
    overlaps = overlaps(mask, :);
    lengths = lengths(mask);
end

function [seqs, overlaps, lengths] = removeNonLetters(seqs, overlaps, lengths)
    fprintf('non letters\n');
    mask = max(seqs, [], 2) <= 4;
    [seqs, overlaps, lengths] = reduceData(mask, seqs, overlaps, lengths);
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

function [seqs, overlaps, lengths] = removeLowHeight(seqs, overlaps, lengths, numerOfTissues, topPercent)
    fprintf('height\n');
    % remove low peaks
    mask = false(length(seqs), 1);
    for i = 1:numerOfTissues
        [vals, ind] = sort(overlaps(:,i), 'descend');
        amountToKeep = round(sum(overlaps(:,i) > 0)*topPercent);
        if vals(1) == vals(end)
            % not an enhancer height peak, either background or gene peak
            amountToKeep = sum(overlaps(:,i) > 0);
        end
        mask(ind(1:amountToKeep)) = true;
    end
    [seqs, overlaps, lengths] = reduceData(mask, seqs, overlaps, lengths);
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
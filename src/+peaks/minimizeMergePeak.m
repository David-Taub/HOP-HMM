% fields of minimizeMergePeak - ['seqs', 'overlaps', 'lengths', 'tissueNames']
function mergedPeaksMin = minimizeMergePeak(mergedPeaks, tissueNames, outFilepath, topPercent)
    L = length(mergedPeaks(1).seq);
    assert(length(mergedPeaks(end).seq) == L);
    % mergedPeaks = removeNonDistal(mergedPeaks);
    [overlaps, lengths] = extractOverlaps(mergedPeaks);
    size(overlaps, 1)
    seqs = trimSeqs(mergedPeaks, L);
    size(overlaps, 1)
    [seqs, overlaps, lengths] = removeLowHeight(seqs, overlaps, lengths, topPercent);
    size(overlaps, 1)
    [seqs, overlaps, lengths] = removeNonLetters(seqs, overlaps, lengths);
    size(overlaps, 1)
    % [seqs, overlaps, lengths] = balanceOverlaps(seqs, overlaps, lengths);
    % outFilepath = '../data/peaks/mergedPeaksMinimized.mat';
    save(outFilepath, '-v7.3', 'seqs', 'overlaps', 'lengths', 'tissueNames');
    fprintf('Saved peaks in %s\n', outFilepath);
    mergedPeaksMin.seqs = seqs;
    mergedPeaksMin.overlaps = overlaps;
    mergedPeaksMin.lengths = lengths;
    mergedPeaksMin.tissueNames = tissueNames;
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
        center = round(length(seq) / 2);
        % seqs(i, :) = nt2int(seq(center-L/2+1:center+L/2));
        start_pos = center - L / 2 +1;
        end_pos = start_pos + L - 1 ;
        seqs(i, :) = seq(start_pos:end_pos);
    end
end


% remove sequences with peaks values that are not inside the top percentage in any tissue
function [seqs, overlaps, lengths] = removeLowHeight(seqs, overlaps, lengths, topPercent)
    fprintf('height\n');
    numerOfTissues = size(overlaps, 2);
    % remove low peaks
    mask = false(length(seqs), 1);
    for i = 1:numerOfTissues
        [vals, ind] = sort(overlaps(:, i), 'descend');
        amountToKeep = round(sum(overlaps(:, i) > 0) * topPercent);
        if vals(1) == vals(end)
            % not an enhancer height peak, either background or gene peak
            amountToKeep = sum(overlaps(:, i) > 0);
        end
        mask(ind(1:amountToKeep)) = true;
    end
    [seqs, overlaps, lengths] = reduceData(mask, seqs, overlaps, lengths);
end


function [overlaps, lengths] = extractOverlaps(mergedPeaks)
    fprintf('overlaps\n');
    numerOfTissues = length(mergedPeaks(1).overlap);
    overlapsFlat = [mergedPeaks.overlap];
    lengths = [mergedPeaks.peakLength]';
    % N x numerOfTissues
    overlaps = reshape(overlapsFlat, [numerOfTissues, length(mergedPeaks)])';
end


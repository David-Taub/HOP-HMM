% fields of minimizeMergePeak - ['seqs', 'overlaps', 'peakLengths', 'tissueNames']
function mergedPeaksMin = minimizeMergePeak(topPercent, doEnhSpecific, withBackground, withGenes,...
                                            seqsPerTissue, L, peakMinL, peakMaxL)
    minimizedMergedFilePath = sprintf('../data/peaks/mergedPeaksMinimized_L%db%dg%dp%des%dspt%dmin%dmax%d.mat', L, ...
                                      withBackground, withGenes, floor(100 * topPercent), doEnhSpecific, seqsPerTissue, ...
                                      peakMinL, peakMaxL);
    if isfile(minimizedMergedFilePath)
        mergedPeaksMin = load(minimizedMergedFilePath);
        return
    end
    [mergedPeaks, tissueNames] = peaks.mergePeakFiles(withBackground, withGenes, true, mergedFilePath, L);
    assert(length(mergedPeaks(end).seq) == L);
    [overlaps, peakLengths, peakPos] = extractOverlaps(mergedPeaks);
    size(overlaps, 1)
    seqs = trimSeqs(mergedPeaks, L);
    size(overlaps, 1)
    [seqs, overlaps, peakLengths, peakPos] = removeLowHeight(seqs, overlaps, peakLengths, peakPos, topPercent);
    size(overlaps, 1)
    [seqs, overlaps, peakLengths, peakPos] = removeNonLetters(seqs, overlaps, peakLengths, peakPos);
    if (peakMaxL > 0) & (peakMinL > 0)
        size(overlaps, 1)
        [seqs, overlaps, peakLengths, peakPos] = limitPeakLength(seqs, overlaps, peakLengths, peakPos, peakMinL, peakMaxL)
    end
    size(overlaps, 1)
    if doEnhSpecific
        [seqs, overlaps, peakLengths, peakPos] = enhancerSpecific(seqs, overlaps, peakLengths, peakPos);
        size(overlaps, 1)
        if seqsPerTissue > 0
            [seqs, overlaps, peakLengths, peakPos] = balanceOverlaps(seqs, overlaps, peakLengths, peakPos, seqsPerTissue);
            size(overlaps, 1)
        end
    end
    % outFilepath = '../data/peaks/mergedPeaksMinimized.mat';
    save(minimizedMergedFilePath, '-v7.3', 'seqs', 'overlaps', 'peakLengths', 'tissueNames', 'peakPos');
    fprintf('Saved peaks in %s\n', outFilepath);
    mergedPeaksMin.seqs = seqs;
    mergedPeaksMin.overlaps = overlaps;
    mergedPeaksMin.peakLengths = peakLengths;
    mergedPeaksMin.peakPos = peakPos;
    mergedPeaksMin.tissueNames = tissueNames;
end


function [seqs, overlaps, peakLengths, peakPos] = enhancerSpecific(seqs, overlaps, peakLengths, peakPos)
    mask = sum(overlaps > 0, 2) == 1;
    [seqs, overlaps, peakLengths, peakPos] = reduceData(mask, seqs, overlaps, peakLengths, peakPos);
end



function [seqs, overlaps, peakLengths, peakPos] = limitPeakLength(seqs, overlaps, peakLengths, peakPos, peakMinL, peakMaxL)
    mask = peakMinL < peakLengths < peakMaxL;
    [seqs, overlaps, peakLengths, peakPos] = reduceData(mask, seqs, overlaps, peakLengths, peakPos);
end


function [seqs, overlaps, peakLengths, peakPos] = balanceOverlaps(seqs, overlaps, peakLengths, peakPos, seqsPerTissue)
    mask = false(size(seqs, 1), 1);
    for i = 1:size(overlaps, 2)
        [vals, inds] = sort(overlaps(:, i), 1, 'descend');
        assert(vals(seqsPerTissue) > 0);
        mask(inds(1:seqsPerTissue)) = true;
    end
    [seqs, overlaps, peakLengths, peakPos] = reduceData(mask, seqs, overlaps, peakLengths, peakPos);
end


function [seqs, overlaps, peakLengths, peakPos] = reduceData(mask, seqs, overlaps, peakLengths, peakPos)
    seqs = seqs(mask, :);
    overlaps = overlaps(mask, :);
    peakLengths = peakLengths(mask);
end


function [seqs, overlaps, peakLengths, peakPos] = removeNonLetters(seqs, overlaps, peakLengths, peakPos)
    fprintf('non letters\n');
    mask = max(seqs, [], 2) <= 4;
    [seqs, overlaps, peakLengths, peakPos] = reduceData(mask, seqs, overlaps, peakLengths, peakPos);
end


function seqs = trimSeqs(mergedPeaks, L)
    fprintf('sequence\n');
    seqsCells = {mergedPeaks.seq};
    seqs = zeros(length(mergedPeaks), L);
    for i = 1:length(seqsCells)
        seq = seqsCells{i};
        center = round(length(seq) / 2);
        % seqs(i, :) = nt2int(seq(center-L/2+1:center+L/2));
        start_peakPos = center - L / 2 +1;
        end_peakPos = start_peakPos + L - 1 ;
        seqs(i, :) = seq(start_peakPos:end_peakPos);
    end
end


% remove sequences with peaks values that are not inside the top percentage in any tissue
function [seqs, overlaps, peakLengths, peakPos] = removeLowHeight(seqs, overlaps, peakLengths, peakPos, topPercent)
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
    [seqs, overlaps, peakLengths, peakPos] = reduceData(mask, seqs, overlaps, peakLengths, peakPos);
end


function [overlaps, peakLengths, peakPos] = extractOverlaps(mergedPeaks)
    fprintf('overlaps\n');
    numerOfTissues = length(mergedPeaks(1).overlap);
    overlapsFlat = [mergedPeaks.overlap];
    peakLengths = [mergedPeaks.peakLength]';
    peakPos = [mergedPeaks.peakpeakPos]' - [mergedPeaks.seqFrom]';
    % N x numerOfTissues
    overlaps = reshape(overlapsFlat, [numerOfTissues, length(mergedPeaks)])';
end

